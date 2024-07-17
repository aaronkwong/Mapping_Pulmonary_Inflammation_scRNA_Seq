suppressPackageStartupMessages(library(optparse))
source("gdpath.R")
source("./R_Scripts/utilities/path_manipulation.R")


option_list <- list(
    make_option(c("-r", "--seurat_object_path"), default="NA",
        help=""),
  #
    make_option(c("-q", "--cluster_resolution"), default="NA",
        help=""),
  #
    make_option(c("-s", "--include_recip"), default="NA",
        help=""),
  #
    make_option(c("-x", "--min_geneset"), default="NA",
        help="The smallest gene sets to be used."),
  #
    make_option(c("-y", "--max_geneset"), default="NA",
        help="The largest gene sets to be used."),
  #
    make_option(c("-t", "--root"), default="NA",
        help=""),  
  #
    make_option(c("-p", "--root_external"), default="NA",
        help=""),
  #
    make_option(c("-c", "--project_name"), default="NA",
        help=""),
  #
    make_option(c("-z", "--logFC_thresh"), default="NA",
        help="")
  )

#load libs and args
opt <- parse_args(OptionParser(option_list=option_list))
options(future.globals.maxSize = 120000 * 1024^2)
library(doParallel)
registerDoParallel(cores=8)
library(Seurat)
source("./R_Scripts/combined_analysis/differential_gene_expression_functions")
source("./R_Scripts/combined_analysis/pathway_analysis_functions.R")

#set vars
seurat_object_path<-opt$seurat_object_path
cluster_resolution<-opt$cluster_resolution
max_geneset<-as.numeric(opt$max_geneset)
min_geneset<-as.numeric(opt$min_geneset)
root<-opt$root
geneset_database_pathways<-c(opt$geneset_database_pathways_1,opt$geneset_database_pathways_2,opt$geneset_database_pathways_3)
geneset_database_tf<-c(opt$geneset_database_tf_1,opt$geneset_database_tf_2,opt$geneset_database_tf_3)
logFC_thresh<-as.numeric(opt$logFC_thresh)

#check whether recipient cells should be included in analysis
if (opt$include_recip=="TRUE"){
  include_recip<-TRUE
}else{
  include_recip<-FALSE
}

#whether files are being written locally or to a cloud synced folder
if (opt$root_external=="TRUE"){
  root_external<-TRUE
}else{
  root_external<-FALSE
}
project_name<-opt$project_name

#create directory where differential gene expression and pathway analysis will be written out
subroot<-paste0(root,"enrichment_analysis/")

#read in the preprocesed object
seurat_object<-readRDS(seurat_object_path)

#the pathway analysis function takes the seurat object and for every cluster in the metadata column called "integrated_clusters" computes differential gene expression between CIT and 2hR timepoints.
#to reproduce the differential gene expression results in the paper of donor cell types between CIT and 2hR we should use the "_lowcluster" metadata
if (cluster_resolution=="low"){
    seurat_object@meta.data[,"integrated_clusters"]<-seurat_object@meta.data[,"integrated_clusters_lowcluster"]
    seurat_object@meta.data[,"integrated_clusters_recipient_origin"]<-seurat_object@meta.data[,"integrated_clusters_recipient_origin_lowcluster"]
    seurat_object@meta.data[,"GSVA_anno_labels"]<-seurat_object@meta.data[,"GSVA_anno_labels_lowcluster"]
    seurat_object@meta.data[,"GSVA_anno_labels_nickname"]<-seurat_object@meta.data[,"GSVA_anno_labels_nickname_lowcluster"]
}else{
    stop("Error. This analysis should only be at low resolution")
}

#should recipient cells be included in the analysis?
if (include_recip==FALSE){
    seurat_object<-seurat_object[,seurat_object@meta.data$recipient_origin==FALSE]
}

#run pathway analysis module which computes differential gene expression and pathway analysis
rank_files<-run_pathway_analysis_combined( 
  seurat_object=seurat_object,
  subroot=subroot,
  root=root,
  root_external=root_external,
  project_name=project_name,
  thresh=thresh
)

saveRDS(rank_files,file=gdpath("rnk_file_names.rds",root=subroot,external=root_external))
