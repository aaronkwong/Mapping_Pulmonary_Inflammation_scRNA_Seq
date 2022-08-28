suppressPackageStartupMessages(library(future))
plan("multiprocess", workers = availableCores())
library(DoubletFinder)

# these functions are modified from:
# https://github.com/chris-mcginnis-ucsf/DoubletFinder

# custom version of paramsweep function where a range of pN and pK is tested. Modified to have a greater range
paramSweep_v3_aw <- function(seu, PCs=1:10, sct = FALSE, num.cores=1) {
  require(Seurat); require(fields);
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.4,by=0.05)

  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]

  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands

  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <- seu@assays$RNA@counts[ , real.cells]
    n.real.cells <- ncol(data)
  }

  if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }

  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep_v3,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,mc.cores=num.cores)
    stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep_v3,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }

  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }

  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)

}

# this function will find and annotate doublets as a new meta data column in the seurat object
# IF parameters for the doublet removal are not specified, a range of values are swept and written
# to help the user select the best parameters
# if the input manually_selected_pN=0 and manually_selected_pK=0 this indicates this is the first run and that just the sweep should be done
# if the manually_selected_pN and manually_selected_pK do not equal 0, then it means they are provdided and the actual doublet removal should be done
# paste0(project_name,"_",sample_name,"_doublet_removal_sweeps.rdata") is the saved sweep results
# paste0(project_name,"_",sample_name,"_doublet_removal.rdata" is the doublet removed processed object
run_find_doublets<-function(seurat_object_processed,pcs,poission_doublet_rate,project_name,sample_name,manually_selected_pN,manually_selected_pK,root,root_external){
	#make pdf figures
	if(manually_selected_pN==0 & manually_selected_pK==0){
    cat(paste0("\npN and pK set to 0. Doing first sweeps as it is the first run.\n"))
		pdf(gdpath(paste0(project_name,"_",sample_name,"_sweeps",".pdf"),root=root,external=root_external))
		#CUSTOM SWEEP USED HERE
		sweep.res.seurat_object_processed <- paramSweep_v3_aw(seurat_object_processed, PCs = 1:pcs, sct = TRUE)
		sweep.stats.seurat_object_processed <- summarizeSweep(sweep.res.seurat_object_processed, GT = FALSE)
		bcmvn_seurat_object_processed <- find.pK(sweep.stats.seurat_object_processed)
		plot(levels(bcmvn_seurat_object_processed$pK)[bcmvn_seurat_object_processed$pK],bcmvn_seurat_object_processed$BCmetric,type="l",col="black",xlab="pK",ylab="BCmvn",main=paste0(sample_name," pK Sweep"))
    save(sweep.res.seurat_object_processed,sweep.stats.seurat_object_processed,bcmvn_seurat_object_processed,file=gdpath(paste0(project_name,"_",sample_name,"_doublet_removal_sweeps.rdata"),root=root,external=root_external))
		dev.off()
    sink(gdpath(paste0(project_name,"_",sample_name,"_sweeps_max",".txt"),root=root,external=root_external))
    print(paste0("max pK ",bcmvn_seurat_object_processed$pK[which.max(bcmvn_seurat_object_processed$BCmetric)], " max BCmetric ",max(bcmvn_seurat_object_processed$BCmetric)))
    sink()
	}else{
    cat(paste0("\npN and pK detected. Continuing with doublet removal.\n"))
		pdf(gdpath(paste0(project_name,"_",sample_name,".pdf"),root=root,external=root_external))
		## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
		homotypic.prop <- modelHomotypic(seurat_object_processed$integrated_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
		nExp_poi <- round(poission_doublet_rate*length(seurat_object_processed$orig.ident))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
		#calculate the number of doublets poisson (homo and hetero) minus the homo doublets
		nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

		## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
		# doublet finder is insensitive to doublets of the same cell type, whereas it should detect doublets of two distinct cell types
    # therefore, the number of expected doublets we should tell doubletFinder to look for is poisson distribution (from 10X technical sheeets)
    # subtracted by the estimated homotypic doublets (homotypic.prop from above). This difference (nExp_poi.adj) are the number of doublets that doubletFinder
    # should find.
    seurat_object_processed <- doubletFinder_v3(seurat_object_processed, PCs = 1:pcs, pN = manually_selected_pN, pK = manually_selected_pK, nExp = nExp_poi.adj, reuse.pANN = FALSE,sct=TRUE)
		#seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = TRUE)
		print(DimPlot(seurat_object_processed , group.by = names(seurat_object_processed@meta.data)[length(names(seurat_object_processed@meta.data))], combine = FALSE,label=T))
		dev.off()
  	save(homotypic.prop,nExp_poi,seurat_object_processed,seurat_object_processed,file=gdpath(paste0(project_name,"_",sample_name,"_doublet_removal.rdata"),root=root,external=root_external))
  	return(seurat_object_processed)
  }
}

