# This script takes an integrated seurat object (and the clustering information) and deconvolves the background RNA levels as described by "Decontamination of ambient RNA in single-cell RNA-seq with DecontX" by Yang et al. Outputs a seurat object with background corrected expression.

library(renv)
renv::activate()

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(celda))
options(future.globals.maxSize = 60000 * 1024^2)
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(Seurat))
library(trqwe)
library(gridExtra)

source("gdpath.R")

plan("multiprocess", workers = availableCores())


option_list <- list(
  make_option(c("-r", "--root"), default="NA",
              help="The directory where the various result directories for each part of the pipeline will be written to."),
  #
  make_option(c("-t", "--root_external"), default="NA",
              help="TRUE or FALSE. If the root is contatined in the cloud synced folder."),  
  #
  make_option(c("-b", "--project_name"), default="NA",
              help="A name for the project.")
  )

opt <- parse_args(OptionParser(option_list=option_list))

#read in run parameters
root<-opt$root
#convert text input to logical
if(opt$root_external=="TRUE"){
    root_external<-TRUE
}else if(opt$root_external=="FALSE"){
    root_external<-FALSE
}else{
  cat("root_external was not defined as true or flase.\n")
}
project_name<-opt$project_name

#read in saved data file 
cat("reading in combine data...\n")
combined_data<-mcreadRDS(gdpath(paste0("initial_integration/combined_data",".rds"),root=root,external=root_external))

cat("creating results folder.\n")
subroot<-paste0(root,"background_removal/")
dir.create(subroot)

#lets look at well known marker genes in our dataset (AT1, AT2, T Cells) before decon, we will also check after
genes_to_examine<-c("SFTPC","AGER","CD3D")
prior_to_decon<-list()
i<-1
for (gene in genes_to_examine){
    prior_to_decon[[i]]<-VlnPlot(combined_data,gene,assay="RNA",slot="counts")
    i<-i+1
}

cat(paste0("running deconvolution with batches",table(combined_data@meta.data$ltx_case,combined_data@meta.data$timepoint)," as batch data.\n"))

#set batch to be each timepoint x case so that every sample is considered its own batch
decontX_data<-decontX(x=combined_data@assays$RNA@counts,z=combined_data@meta.data$integrated_clusters,batch=paste0(combined_data@meta.data$ltx_case,"_",combined_data@meta.data$timepoint))


#we should make some plots showing the identified contamination
pdf(gdpath("decontamination_plots.pdf",root=subroot,external=root_external))
combined_data <- AddMetaData(
    object = combined_data,
    metadata = decontX_data$contamination,
    col.name = 'decontX_contamination'
    )
FeaturePlot(combined_data,features="decontX_contamination")
dev.off()

#confirm that deconvolution has altered the matrix
if(!(identical(decontX_data[["decontXcounts"]],combined_data@assays$RNA@counts))){
    cat("decontaminated counts are indeed different from original.\n")
}else{
    cat("error? decontaminated counts are the same as original??")
    quit()
}

#the @count and @data slot in assay RNA should be identical so well replace both
if(identical(combined_data@assays$RNA@counts,combined_data@assays$RNA@data)){
  print("replacing @assays$RNA@counts and @assays$RNA@data")
  combined_data@assays$RNA@counts<-decontX_data[["decontXcounts"]]
  combined_data@assays$RNA@data<-decontX_data[["decontXcounts"]]
}else{
  cat("original RNA slots counts and data are not identical. Fix this before explacing with decontX values")
}

#saving the decontX_data file
mcsaveRDS(decontX_data,file=gdpath("decontX_data.rds",root=subroot,external=root_external))
mcsaveRDS(combined_data,file=gdpath("combined_data_decontaminated.rds",root=subroot,external=root_external))

#save some stats on the final product
sink(gdpath("combined_data_summary.txt",root=subroot,external=root_external))
combined_data
sink()


#lets look at expression after background removal to see if markers become more strongly expressed
pdf(gdpath("pre_post_decont_cluster_expression.pdf",root=subroot,external=root_external))
for (i in 1: length(genes_to_examine)){
    grid.arrange(grobs=list(prior_to_decon[[i]]+ ggtitle(paste0("pre_decontamination_",genes_to_examine[i])),after_decon[[i]] + ggtitle(paste0("post_decontamination_",genes_to_examine[i]))),ncol=2)
}
dev.off()

cat("Background removal complete")
quit()