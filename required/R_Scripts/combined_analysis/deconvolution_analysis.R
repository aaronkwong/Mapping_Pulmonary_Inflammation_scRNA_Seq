suppressPackageStartupMessages(library(renv))
renv::activate()
suppressPackageStartupMessages(library(optparse))
library(Seurat)
library(trqwe)
source("gdpath.R")
source("./R_Scripts/combined_analysis/pathway_analysis_functions.R")
source("./R_Scripts/combined_analysis/differential_gene_expression_functions.R")


option_list <- list(
    make_option(c("-p", "--seurat_object_path"), default="NA",
		help=""),
  #
    make_option(c("-x", "--cluster_resolution"), default="NA",
        help=""),
  #
    make_option(c("-y", "--include_recip"), default="NA",
        help=""),
    make_option(c("-a", "--case_list"), default="NA",
              help=""),
    make_option(c("-e", "--clusters_with_recip_cells"), default="NA",
              help="cluster numbers containing recip cells seperated by commas"),
    make_option(c("-v", "--root"), default="NA",
              help=""),
    make_option(c("-q", "--root_external"), default="NA",
              help=""),
	make_option(c("-r", "--analysis_name"), default="NA",
              help="")

  )

opt <- parse_args(OptionParser(option_list=option_list))
cluster_resolution<-opt$cluster_resolution
analysis_name<-opt$analysis_name


#extract arguments
seurat_object_path<-opt$seurat_object_path
root<-opt$root
case_list<-unlist(strsplit(opt$case_list,split=","))
geneset_database_pathways<-opt$geneset_database_pathways
geneset_database_tf<-opt$geneset_database_tf
recip_containing_cluser<-as.numeric(unlist(strsplit(opt$clusters_with_recip_cells,split=",")))

case_name<-paste0(case_list,collapse="_")

#check whether recipient cells should be included in analysis
if (opt$include_recip=="TRUE"){
  include_recip<-TRUE
}else{
  include_recip<-FALSE
}


# parameters for writing output
if(opt$root_external=="TRUE"){
  root_external<-TRUE
}else{
  root_external<-FALSE
  cat("internal root detected.\n")
}


#make output directory results directories
results_dir<-paste0(root,"deconvolution_analysis_",analysis_name,"/")
dir.create(results_dir)
cat("Creating output results directories...\n")
subroot<-paste0(results_dir,"pathway_analysis/")
case_subroot<-paste0(subroot,paste0(case_list,collapse="_"),"/")
suppressWarnings(dir.create(subroot))
suppressWarnings(dir.create(case_subroot))

cat("reading in seurat object...\n")
seurat_object<-mcreadRDS(seurat_object_path)

#all functions below this line read the clustering data from the metadata column "integrated_clusters"
#Populate the "integrated_clusters" with either the "_lowcluster" (low resolution) or "_subcluster" (high resolution)
if (cluster_resolution=="low"){
    seurat_object@meta.data[,"integrated_clusters"]<-seurat_object@meta.data[,"integrated_clusters_lowcluster"]
    seurat_object@meta.data[,"integrated_clusters_recipient_origin"]<-seurat_object@meta.data[,"integrated_clusters_recipient_origin_lowcluster"]
    seurat_object@meta.data[,"GSVA_anno_labels"]<-seurat_object@meta.data[,"GSVA_anno_labels_lowcluster"]
    seurat_object@meta.data[,"GSVA_anno_labels_nickname"]<-seurat_object@meta.data[,"GSVA_anno_labels_nickname_lowcluster"]
}else if (cluster_resolution=="high"){
    seurat_object@meta.data[,"integrated_clusters"]<-seurat_object@meta.data[,"integrated_clusters_subcluster"]
    seurat_object@meta.data[,"integrated_clusters_recipient_origin"]<-seurat_object@meta.data[,"integrated_clusters_recipient_origin_subcluster"]
    seurat_object@meta.data[,"GSVA_anno_labels"]<-seurat_object@meta.data[,"GSVA_anno_labels_subcluster"]
    seurat_object@meta.data[,"GSVA_anno_labels_nickname"]<-seurat_object@meta.data[,"GSVA_anno_labels_nickname_subcluster"]
}else{
	stop("error cluster_resolution must be either low or high")
}

#should recipient cells be included in the analysis? if False remove all recipient cells
if (include_recip==FALSE){
    seurat_object<-seurat_object[,seurat_object@meta.data$recipient_origin==FALSE]
}

#pathway analysis
rank_files<-c()
rank_files_names<-c()
de_result<-list()
de_result_names<-c()
x<-1
for (cell_type in recip_containing_cluser){
	skipped_de<-0
	recip_cells_post<-rownames(seurat_object@meta.data)[(seurat_object@meta.data$integrated_clusters_recipient_origin==paste0(cell_type,"_","TRUE")) & (seurat_object@meta.data$timepoint=="post")]
	donor_cells_post<-rownames(seurat_object@meta.data)[(seurat_object@meta.data$integrated_clusters_recipient_origin==paste0(cell_type,"_","FALSE")) & (seurat_object@meta.data$timepoint=="post")]
	donor_cells_pre<-rownames(seurat_object@meta.data)[(seurat_object@meta.data$integrated_clusters_recipient_origin==paste0(cell_type,"_","FALSE")) & (seurat_object@meta.data$timepoint=="pre")]
	if(file.exists(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_pre.rds"),root=case_subroot,external=root_external)) & file.exists(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_post.rds"),root=case_subroot,external=root_external))){
		cat(paste0("found existing save files. Loading them...\n",cell_type,"\n"))
		compare_recip_to_donor_pre<-mcreadRDS(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_pre.rds"),root=case_subroot,external=root_external))
		compare_recip_to_donor_post<-mcreadRDS(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_post.rds"),root=case_subroot,external=root_external))
	}else{
		cat("running differential gene expression...\n")
		if(length(recip_cells_post)>20){
			#conduct differential expression
			compare_recip_to_donor_pre<-run_nebula_raw_counts_deconvolution(seurat_object=seurat_object,ident.1=donor_cells_pre, ident.2=recip_cells_post,metadata=metadata)
			compare_recip_to_donor_post<-run_nebula_raw_counts_deconvolution(seurat_object=seurat_object,ident.1=donor_cells_post, ident.2=recip_cells_post,metadata=metadata)
			mcsaveRDS(compare_recip_to_donor_pre,gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_pre.rds"),root=case_subroot,external=root_external))
			mcsaveRDS(compare_recip_to_donor_post,gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_post.rds"),root=case_subroot,external=root_external))
		}else{
			cat(paste0("\n\n###################################################################\n"))
			cat(paste0(cell_type," did not have at least 50 cells in the recipient group.\n"))
			cat(paste0("###################################################################\n"))
			skipped_de<-1
		}
	}
	if (skipped_de!=1){
		de_result[[x]]<-compare_recip_to_donor_pre
		de_result[[x+1]]<-compare_recip_to_donor_post
		x<-x+2
		de_result_names<-c(de_result_names,paste0(cell_type,c("_recip_to_donor_pre","_recip_to_donor_post")))
		#make histogram plots
		cat("making histogram plots\n")
		pdf(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","pvalue_distribution.pdf"),root=case_subroot,external=root_external))
		print(hist(compare_recip_to_donor_pre[,"p_val"],breaks=20,main=paste0(cell_type,": recipient_post vs donor_pre")))
		print(hist(compare_recip_to_donor_post[,"p_val"],breaks=20,main=paste0(cell_type,": recipient_post vs donor_post")))
		dev.off()

		#write tab files
		write.table(compare_recip_to_donor_pre,file=gdpath(paste0(cell_type,"_recipient_post_vs_donor_pre_DE",".tab"),root=case_subroot,external=root_external),sep="\t",row.names=TRUE,col.names=NA)
		write.table(compare_recip_to_donor_post,file=gdpath(paste0(cell_type,"_recipient_post_vs_donor_post_DE",".tab"),root=case_subroot,external=root_external),sep="\t",row.names=TRUE,col.names=NA)	
	}
}

names(de_result)<-de_result_names
mcsaveRDS(de_result,gdpath(paste0("all_comparisons_recip_to_donor_post.rds"),root=case_subroot,external=root_external))
