# This script takes the results from the background removal and doublet removal and integrates the data. This is done by taking the seurat_object containing the background corrected gene expression, and then further removing any cells which were called as doublets (from output of doubletFinder). Note: background corrected expression is rounded to the nearest integer. Background removal is done on a per-sample basis. Therefore, all samples are re-integrated as per "Comprehensive Integration of Single-Cell Data" by Stuart et al. This script outputs a single integrated seurat object and some summary plots.

library(renv)
renv::activate()

suppressPackageStartupMessages(library(optparse))
library(Seurat)
library(trqwe)
source("gdpath.R")
suppressPackageStartupMessages(library(future))
options(future.globals.maxSize = 1200000 * 1024^2)
plan("multiprocess", workers = 16)

option_list <- list(
  make_option(c("-r", "--root"), default="NA",
              help="he directory where the various result directories for each part of the pipeline will be written to."),
  #
  make_option(c("-t", "--root_external"), default="NA",
              help="TRUE or FALSE. If the root is contatined in the cloud synced folder."),  
  #
  make_option(c("-a", "--pc_use"), default="NA",
              help="The number of principal components to use when generating UMAP."),
  #
  make_option(c("-b", "--project_name"), default="NA",
              help="A name for the project."),
  #
  make_option(c("-d", "--project_data_path"), default="NA",
              help="The directory containing a .tab file containing the information for each sample."),
  #
  make_option(c("-c", "--cluster_resolution"), default="NA",
              help="clustering resolution to run the analysis at."),
  #
  make_option(c("-f", "--remove_doublets"), default="TRUE",
              help="TRUE or FALSE. Should cells annotated as doublets be removed?")
  )

opt <- parse_args(OptionParser(option_list=option_list))

root<-opt$root
project_data_path<-opt$project_data_path
project_data<-read.csv(project_data_path,header=TRUE)
project_working<-project_data[project_data$directory!="",]
#convert text input to logical for root_external
if(opt$root_external=="TRUE"){
    root_external<-TRUE
}else if(opt$root_external=="FALSE"){
    root_external<-FALSE
}else{
  cat("root_external was not defined as true or flase.\n")
}
#convert text input to logical for remove doublets
if(opt$remove_doublets=="TRUE"){
    remove_doublets<-TRUE
}else if(opt$remove_doublets=="FALSE"){
    remove_doublets<-FALSE
}else{
  cat("remove_doublets was not defined as true or flase.\n")
}
pc_use<-as.numeric(opt$pc_use)
project_name<-opt$project_name
cluster_resolution=as.numeric(opt$cluster_resolution)

source(gdpath("Masters/project_6_scRNA_seq/scRNA_seq_Cell_Tracker/combined_analysis/integration_functions.R"))
source(gdpath("Masters/project_6_scRNA_seq/scRNA_seq_Cell_Tracker/annotate_cells/annotate_cells.R"))
cat("running rPCA anchored integration.\n")

print(project_data)

subroot<-paste0(root,"post_QC_processing/")
dir.create(subroot)

#load doublet and decontX data
cat("reading in doublet and background removed data.\n")
combined_data_doublet_finder<-mcreadRDS(gdpath("combined_data_doublet_finder.rds",root=paste0(root,"doublet_removal/"),external=root_external))
combined_data_decontaminated<-mcreadRDS(gdpath("combined_data_decontaminated.rds",root=paste0(root,"background_removal/"),external=root_external))

#now lets clean up the doublet data
# if we look at the doublet data meta data we get many classification columns, becase we merged the indiviual items.
# if we take all the columns in met data with the "DF.classification" in them then we have a table
# from that table we need to take the non NA value which would be either "Singlet" or "Doublet"
# [15] "pANN_0.25_0.11_1211"                "DF.classifications_0.25_0.11_1211"
# [17] "pANN_0.25_0.01_33"                  "DF.classifications_0.25_0.01_33"
# [19] "pANN_0.25_0.04_70"                  "DF.classifications_0.25_0.04_70"
# [21] "pANN_0.25_0.001_1052"               "DF.classifications_0.25_0.001_1052"

# we need to append the metadata of the doublet finder analysis to the decontaminated seurat object
if(all(rownames(combined_data_doublet_finder@meta.data)==rownames(combined_data_decontaminated@meta.data))){
    cat("barcodes of cells from decontX and doubletFinder match. Appending classifications from DoubletFinder to DecontX object\n")
    #now lets append this doublet data to the decontamination seurat object
    combined_data <- AddMetaData(
        object = combined_data_decontaminated,
        metadata = combined_data_doublet_finder$DoubletFinderClassification,
        col.name = 'doublet_finder'
    )
    combined_data <- AddMetaData(
        object = combined_data,
        metadata = combined_data_doublet_finder$DoubletFinderScores,
        col.name = 'doublet_finder_scores'
    )
}else{
    cat("barcodes from the outputs of doubletfinder and DecontX don't match. Please fix.\n")
    quit()
}



#now lets remove the decon and doublet remove objects to save memory
rm(combined_data_doublet_finder)
rm(combined_data_decontaminated)
gc()

if(remove_doublets){
    # now that we have the object with the decontaminated matrix and doublet calls# lets remove the cells called as doublets
    cat(paste0("\nThere were ",nrow(combined_data@meta.data)," cells prior to doublet removal\n"))
    combined_data<-combined_data[,combined_data@meta.data$doublet_finder=="Singlet"]
    cat(paste0("There were ",nrow(combined_data@meta.data)," cells after to doublet removal\n"))
}else{
    cat("remove doublets is FALSE. Doublets scores and classification will be added to meta data, but no cells will be removed.\n")
}

#set deauflt assay back to RNA and set defualt assay to RNA
DefaultAssay(combined_data)<-"RNA"
combined_data[["SCT"]]<-NULL

unique_gene_cutoff<-200
#decontX gives non integer values. Round these to integers for other analysis that expect integer values only.
combined_data@assays$RNA@counts<-round(combined_data@assays$RNA@counts)
#since decontX changed the raw count matrix (also rounded), the nFeature and nCount in seurat is no longer accurate. Recalculate
combined_data$nFeature_RNA<-colSums(combined_data@assays$RNA@counts>0)
combined_data$nCount_RNA<-colSums(combined_data@assays$RNA@counts)
cat(paste0(sum(combined_data$nFeature_RNA<unique_gene_cutoff)," cell removed as they had less than ",unique_gene_cutoff," unqiue genes.\n"))
#only keep cells that have at least 200 unique genes
combined_data<-combined_data[,combined_data$nFeature_RNA>200]

#now lets remove cells called as doublets
cat("\n##############################\n")
cat("removing cells from gavins analysis called as doublets\n")
cat("##############################\n\n")
cat("seurat object has ",nrow(combined_data@meta.data),"cells prior to removal of doublets identified by deconvolution.\n")
gavin_doublets<-c()
for (i in 1:nrow(project_working)){
    if (project_working$deconvolution_path[i]!="no_decon"){
        #extract barcodes from gavin with genotype 2 (mixed, likely doublets)
        doublet_barcodes<-extract_recip_barcodes_from_gavin_table(path=project_working$deconvolution_path[i],label_in_label_column=2) # doublets are alwyas class 2
        cat(paste0("reading deconvolution data from: ",project_working$deconvolution_path[i],".",length(doublet_barcodes)," barcodes found in gavin's data.\n"))
        doublet_barcodes<-paste0(project_working$case[i],project_working$timepoint[i],"_",doublet_barcodes)
        gavin_doublets<-c(gavin_doublets,doublet_barcodes)
    }
}
combined_data<-remove_cells(seurat_object=combined_data,barcodes_of_cells_to_be_removed=gavin_doublets)
cat("seurat object has ",nrow(combined_data@meta.data),"cells after removal of doublets identified by deconvolution.\n")


#for the anchored integration we need to seperate the combiend object back into its individual components. 
all_samples<-list()
i<-1
for (sample_name in project_working$sample_name){
    all_samples[[i]]<-combined_data[,combined_data@meta.data$sample_name==sample_name]
    i<-i+1
}



###################################################################################################
# Start of the Anchored Integrated  normalization technique
###################################################################################################
cat("Starting Anchored Integration of two timepoints....","\n")
#combine two datasets into a single object
# combined_data_list <- c(pre.evlp,post.evlp)

#normalize each dataset
combined_data_list<-list()
for (i in 1:length(all_samples)) {
  # combined_data_list[[i]] <- SCTransform(all_samples[[i]], verbose = T,conserve.memory=TRUE,method = 'glmGamPoi',residual_type="pearson",min_cells=10)
  # combined_data_list[[i]] <- SCTransform(all_samples[[i]], verbose = T,conserve.memory=TRUE,min_cells=10)
  combined_data_list[[i]]  <- SCTransform(all_samples[[i]], verbose = TRUE,conserve.memory=TRUE,min_cells=10)
}

#remove all_samples to free up ram as we no longer need it
rm(all_samples)
gc()

#REPLACED WITH RPCA VERSION
# #calculate pearson features
# combined_data.features <- SelectIntegrationFeatures(object.list = combined_data_list, nfeatures = 3000)
# combined_data_list <- PrepSCTIntegration(object.list = combined_data_list, anchor.features = combined_data.features, verbose = FALSE)

# #calculate anchors
# combined_data.anchors <- FindIntegrationAnchors(object.list = combined_data_list, normalization.method = "SCT", 
#   anchor.features = combined_data.features, verbose = TRUE)
# combined_data.integrated <- IntegrateData(anchorset = combined_data.anchors, normalization.method = "SCT", 
#   verbose = TRUE)
# combined_data.integrated<-run_standard_integration(combined_data_list=combined_data_list,reference_index=c(3))

reference_index=c(3)
combined_data.features <- SelectIntegrationFeatures(object.list = combined_data_list, nfeatures = 2000)
#un PCA for rpca
combined_data_list<-lapply(X = combined_data_list, FUN = function(x) {
    # x <- ScaleData(x, features = features, verbose = FALSE) #can remove RunPCA pulls from scale.data slot of SCT
    x <- RunPCA(x, features = combined_data.features, verbose = FALSE)
})
combined_data_list <- PrepSCTIntegration(object.list = combined_data_list, anchor.features = combined_data.features, verbose = FALSE)
combined_data.anchors <- FindIntegrationAnchors(object.list = combined_data_list,reference=reference_index,anchor.features=combined_data.features, normalization.method = "SCT", verbose = TRUE,reduction="rpca")
rm(combined_data_list)
gc()
combined_data.integrated <- IntegrateData(anchorset = combined_data.anchors, normalization.method = "SCT", 
  verbose = TRUE)


#display umap
combined_data.integrated <- RunPCA(combined_data.integrated, verbose = FALSE)
combined_data.integrated <- RunUMAP(combined_data.integrated, dims = 1:pc_use)
combined_data.integrated <- FindNeighbors(combined_data.integrated, dims = 1:pc_use)
combined_data.integrated <- FindClusters(combined_data.integrated, resolution = cluster_resolution)

#overwrite seurat clusters in the integrated object. This will make seurat_clusters and integrated_clusters identical identical. Useful when plotting.
combined_data.integrated<-AddMetaData(object=combined_data.integrated,metadata=combined_data.integrated@meta.data$seurat_clusters,col.name="integrated_clusters")

#this saves the anchored object with the default assay as anchored. This should only be used cluster label annotation.
mcsaveRDS(combined_data.integrated,file=gdpath(paste0(project_name,"_ANCHORED_integration",".rds"),root=subroot,external=root_external))

cat("Anchored Integration of two timepoints complete.","\n")
###################################################################################################
# End of the Anchored Integrated  normalization technique
###################################################################################################


#read it back now
combined_data<-mcreadRDS(file=gdpath(paste0("combined_data_no_anchored_data",".rds"),root=subroot,external=root_external))

#now we can add the integrated clusters labels on top of the per cell normalization expression data
#first lets make surer the names of the cells in the per cell normalization and the integrated clustering match
if (all(rownames(combined_data@meta.data) == rownames(combined_data.integrated@meta.data))){
    cat(paste0("Transferring labels from integrated clustering onto per-cell normalization\n"))
    combined_data<-AddMetaData(object=combined_data,metadata=combined_data.integrated@meta.data$seurat_clusters,col.name="integrated_clusters")
    #set the integrated clustering as the Idents
    Idents(combined_data)<-combined_data@meta.data$integrated_clusters
}else{
    cat(paste0("Error. Barcodes of integrated and per cell normalization did not match.\n"))
}

#this file contains the per cell normalization (also as defualt assay) and with an additional metadata column "integrated clusters" which has been derived from the integrated clustering
mcsaveRDS(combined_data,file=gdpath(paste0(project_name,"_combined_data",".rds"),root=subroot,external=root_external))


#create some plots to see how early clustering looks for integrated data
cat("Creating pdf figures...\n")
pdf(gdpath("sc_integrated.pdf",root=subroot,external=root_external),height=10,width=17)
plots <- DimPlot(combined_data.integrated, group.by = c("seurat_clusters","timepoint"), combine = FALSE,label=T,raster=FALSE)
print(CombinePlots(plots))
plots <- DimPlot(combined_data.integrated, group.by = c("ltx_case"), combine = FALSE,label=T,raster=FALSE)
print(CombinePlots(plots))
plots <- DimPlot(combined_data.integrated, group.by = c("seurat_clusters"), combine = FALSE,label=T,raster=FALSE)
print(CombinePlots(plots))
dev.off()
