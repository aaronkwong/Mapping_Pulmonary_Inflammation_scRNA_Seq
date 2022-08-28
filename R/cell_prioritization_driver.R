#This script runs the augur cell prioritization analysis on donor cells comparing profiles pre and post transplant to identify the cells most changed. Method is as described in "Cell type prioritization in single-cell data" by Skinnider et al. 

suppressPackageStartupMessages(library(optparse))
source("gdpath.R")
source("augur_functions.R")



option_list <- list(
    make_option(c("-r", "--seurat_object_path"), default="NA",
        help="path to the seurat_object"), 
    make_option(c("-t", "--root"), default="NA",
        help="The directory where the various result directories for each part of the pipeline will be written to."),  
    make_option(c("-p", "--root_external"), default="NA",
        help="TRUE or FALSE. If the root is contatined in the cloud synced folder.")
    )

opt <- parse_args(OptionParser(option_list=option_list))

#get run parameters
seurat_object_path<-opt$seurat_object_path
root<-opt$root
if (opt$root_external=="TRUE"){
  root_external<-TRUE
}else{
  root_external<-FALSE
}

cat("sourcing utility script...\n")


#create results directory
cat("creating results folder...\n")
subroot<-paste0(root,"cell_prioritization/")
dir.create(subroot)

if(dir.exists(subroot)){
}else{
    cat(paste0("Error could not create results directory: ",subroot,"\n"))
}

cat("reading in seurat object...\n")
seurat_object<-readRDS(seurat_object_path)

#run augur analysis
cat("running augur analysis...\n")
run_augur_analysis(seurat_object=seurat_object,root=subroot,root_external=root_external)