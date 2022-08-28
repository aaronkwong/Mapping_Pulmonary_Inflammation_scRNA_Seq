#This script compares a "test" seurat object that needs annotation for its clusters against a known reference dataset (input as another seurat_object) using the package singleR. In this case we will use the Travaglini et al. Lung single cell atlas published in Nature. Script outputs a pdf file showing the top 10 most likely cluster annotations for each cluster in the test set based off of markers from the reference. 

suppressPackageStartupMessages(library(renv))
renv::activate()
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(trqwe))
suppressPackageStartupMessages(library(BiocParallel))
bpparam<-MulticoreParam(16)

option_list <- list(
    make_option(c("-f", "--root"), default="NA",
        help="The directory where the various result directories for each part of the pipeline will be written to."),  
    make_option(c("-g", "--root_external"), default="NA",
        help="TRUE or FALSE. If the root is contatined in the cloud synced folder."),  
    make_option(c("-i", "--travag_reference_path"), default="NA",
        help="path to rds file"),  
    make_option(c("-j", "--test_cells_path"), default="NA",
        help="path to rds file")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# source custom functions
source("gdpath.R")

# root directory
root<-opt$root
# path to molecular lung cell atlas sample data
# path to rds file containing Suerat object with molecular lung cell atlas 
# these cells will be used as reference for annotation
travag_reference_path<-opt$travag_reference_path
#load in the seurat object containing the cells to be annotated against the reference
test_cells_path<-opt$test_cells_path

#create directory to output files
subroot<-paste0(root,"singleR_anno/")
dir.create(subroot)
#get root external
if (opt$root_external=="TRUE"){
  root_external<-TRUE
}else{
  root_external<-FALSE
}

#read in the reference data
cat("reading in reference data...\n")
travag_cells<-mcreadRDS(travag_reference_path)
ref<-travag_cells@assays$SCT@data
ref_anno<-travag_cells@meta.data$free_annotation

#read in the test dataset, or the dataset which you would like to annotate through comparison to the reference
cat("readining test data...\n")
test_cells<-mcreadRDS(test_cells_path)
test<-test_cells@assays$SCT@data

#run singleR
cat("running singleR analysis..\n")
result<-SingleR(
    test=test,
    ref=ref,
    labels=ref_anno,
    BPPARAM=bpparam
)

#make sure cell barcodes of seurat object passed to singleR are still the same after processing
if (!all(rownames(result$scores)==rownames(test_cells@meta.data))){
    cat(paste0("Error. Original object passed to singleR has out of order, or missing cells in the output."))
    quit()
}

# for each cell type print the scores for the top 10 probable cell types annotations and 
# corresponding scores to manual review
pdf(gdpath("cluster_anno.pdf",root=subroot,external=root_external),width=14)
#loop through each cell type and do annotation
for (i in 0:(length(unique(test_cells@meta.data$integrated_clusters))-1)){
	cat("processing cell type: ",i,"\n")
    # only annotate cells using cells from the pre transplant timepoint as this is most close to atlas cells
	if(sum(test_cells@meta.data$integrated_clusters==i & test_cells@meta.data$timepoint=="pre")>0){
		h<-result$scores[test_cells@meta.data$integrated_clusters==i & test_cells@meta.data$timepoint=="pre",]
		cell_type_scores<-colMeans(h)
		cell_type_scores<-cell_type_scores[order(cell_type_scores,decreasing=TRUE)]
		# print(cell_type_scores)
		par(mar=c(2,14,2,2))
		barplot(cell_type_scores[1:10],main=paste0("cluster ",i),las=1,horiz=TRUE,col.axis=1)
	}
}
dev.off()