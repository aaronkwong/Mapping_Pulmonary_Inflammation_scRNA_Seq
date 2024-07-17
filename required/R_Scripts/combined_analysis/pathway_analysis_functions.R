library(Seurat)
library(future)
plan("multiprocess", workers = 12)

# thresh is a number between 0 and 1 defining the minimum proportion of cells that must be expressed in either ident.1 or ident.2 to be kept for analysis
run_pathway_analysis_combined<-function(seurat_object,subroot,root,root_external,project_name,thresh=0,logFC_thresh=0,de_tool="nebula"){
	if( file.exists(gdpath(paste0(project_name,"_all_cluster_DE.rds"),root=subroot,external=root_external))){
		df.cluster<-readRDS(file=gdpath(paste0(project_name,"_all_cluster_DE.rds"),root=subroot,external=root_external))
		rank_files<-readRDS(file=gdpath(paste0(project_name,"DE_rank_files.rds"),root=subroot,external=root_external))
	}else{

		cat("creating pathway_analysis directory output...\n")
		dir.create(subroot)
		cat("starting pathway analysis...\n")
		#create the compound label so that cells of the same cluster but are from different timepoints are labelled accordingly
		compound_labels<-paste0(seurat_object$timepoint,seurat_object$integrated_clusters)
		seurat_object<-AddMetaData(object=seurat_object,metadata=compound_labels,col.name="compound")
		Idents(seurat_object)<-seurat_object$compound
		#lets do a check to make sure that there are enough cells pre and post in each cluster
		if (all(table(seurat_object$timepoint,seurat_object$integrated_clusters)<5)){
			cat("\nError, some clusters have fewer than 5 cells at at least one timepoint: \n")
			print(table(seurat_object$timepoint,seurat_object$integrated_clusters))
			quit()
		}
		#compute differential gene expression for each cluster
		#create vector to store rank file locations
		cat("Calculating cluster specific differential gene expression and making p-value histogram plots...\n")
		pdf(gdpath(paste0("cluster_pathway_analysis_",project_name,".pdf"),root=subroot,external=root_external))
		rank_files<-c()
		df.cluster<-list()
		# df.cluster_name<-c()
		cat(paste0("There are ",length(levels(seurat_object@meta.data$integrated_clusters))," cell types detected.\n"))
		for (i in 0:(length(levels(seurat_object@meta.data$integrated_clusters))-1)) {
		# for (i in 0) {
			#make sure there are enough pre and post cells
			num_pre<-sum(seurat_object@meta.data$compound==paste0("pre",i))
			num_post<-sum(seurat_object@meta.data$compound==paste0("post",i))
			if(num_pre>5 & num_post>5){#if there are enough cells for DE
				if(de_tool=="BPSC"){
					#logfc limit does not benefit by change to -Inf due 
					test.df <- run_cell_type_de(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i,pseudocount.use=1)
					
					#1 must be added since counter starts from 0
					
					de_res<-data.frame(test.df$PVAL,test.df$logFC,p.adjust(test.df$PVAL,method="BH"),test.df$ident.1.pct,test.df$ident.2.pct)
					colnames(de_res)<-c("","avg_logFC","p_val_adj","ident.1.pct","ident.2.pct")
					rownames(de_res)<-names(test.df$PVAL)
					de_res<-de_res[!is.na(de_res$p_val),]
					de_res<-de_res[order(de_res$p_val_adj),]
					de_res<-de_res[(de_res$ident.1.pct>thresh | de_res$ident.2.pct>thresh),]
				}else if(de_tool=="voomlimma"){ #limma voom method
					de_res<-run_cell_type_de_limma_voom(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i)
				}else if(de_tool=="seurat_LR"){ #seurat LR
					de_res<-run_cell_type_de_seurat_LR(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i,thresh=thresh,logFC_thresh=logFC_thresh)
				}else if(de_tool=="seurat_negbinom"){ #seurat negbinom
					de_res<-run_cell_type_de_seurat_negbinom(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i,thresh=thresh,logFC_thresh=logFC_thresh)
				}else if(de_tool=="nebula"){ #nebula
					de_res<-run_nebula_raw_counts(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i)
				}
				df.cluster[[i+1]]<-de_res

				#make historgram
				print(hist(de_res$p_val,breaks=20,main=paste0("cluster_",i,"_min_express_cutoff_",thresh)))
				rank_files<-c(rank_files,gdpath(paste0(project_name,"_cluster_",i,".rnk"),root=subroot,external=root_external))
			}else{
				df.cluster[[i+1]]<-"insufficient_cells"
				rank_files<-c(rank_files,"insufficient_cells")
			}

			# df.cluster_name<-c(df.cluster_name,paste0("cluster_",i))
			print(paste0(i," run complete"))
		}
		dev.off()
		names(rank_files)<-paste0(project_name,"_cluster_",0:as.numeric(length(rank_files)-1))
		names(df.cluster)<-paste0("cluster_",0:as.numeric(length(rank_files)-1))
		#df.cluster is a list object where each list corresponds to the differnetially expressed genes for each cluster
		#rank.files is a list of the names of the rank files generated
		saveRDS(df.cluster,file=gdpath(paste0(project_name,"_all_cluster_DE.rds"),root=subroot,external=root_external))
		saveRDS(rank_files,file=gdpath(paste0(project_name,"DE_rank_files.rds"),root=subroot,external=root_external))

		#save cluster DE
		cat("Saving cluster differential expression....\n")
		for (i in 1:length(df.cluster)){
			write.table(df.cluster[[i]],file=gdpath(paste0(project_name,"_cluster_",i-1,"_DE.tab"),root=subroot,external=root_external),sep="\t",row.names=TRUE,col.names=NA)
		}
		cat("Pathway analysis complete....\n")
	}

	#cycle through cell numbers and grab the DE tables. Then cutoff by thresh (thresh is the minimum proportion of cells in either ident.1 or ident2 that must be met for a gene to be kept)
	for (i in 0:(length(levels(seurat_object@meta.data$integrated_clusters))-1)) {
		if(df.cluster[[i+1]]!="insufficient_cells"){
			subset<-df.cluster[[i+1]][(df.cluster[[i+1]][,"ident.1.pct"]>thresh | df.cluster[[i+1]][,"ident.2.pct"]>thresh),]
			cat("using make_rnk_simple_by_fold_change\n")
			make_rnk_simple(subset[,c("p_val","avg_logFC")],gdpath(paste0(project_name,"_cluster_",i,".rnk"),root=subroot,external=root_external),clean=F)
		}else{
			cat(paste0("cluster ",i," had insufficient cells.\n"))
		}
	}

	return(rank_files)
}


