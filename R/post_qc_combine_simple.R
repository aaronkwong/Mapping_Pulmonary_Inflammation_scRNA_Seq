# This script takes the results from the background removal and doublet removal and integrates the data. This is done by taking the seurat_object containing the background corrected gene expression, and then further removing any cells which were called as doublets (from output of doubletFinder). Note: background corrected expression is rounded to the nearest integer. Background removal is done on a per-sample basis. Therefore, all samples are combined (no batch correction). A saved seruat object and summary plots are generated.  

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
              help="The directory where the various result directories for each part of the pipeline will be written to."),
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

source("integration_functions.R")
source("annotate_cells_functions.R")
cat("running rPCA anchored integration.\n")


unique_gene_cutoff<-200

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
    cat(paste0("There were ",nrow(combined_data@meta.data)," cells after doublet removal\n"))
}else{
    cat("remove doublets is FALSE. Doublets scores and classification will be added to meta data, but no cells will be removed.\n")
}

#set deauflt assay back to RNA and set defualt assay to RNA
DefaultAssay(combined_data)<-"RNA"
combined_data[["SCT"]]<-NULL

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
# Start of the simple per-cell normalization technique
###################################################################################################
cat("Starting simple per-cellnormalization...\n")

#https://github.com/ChristophH/sctransform/blob/master/supplement/batch_correction.Rmd
#if regular SCTransform throws error use the glmGamPoi version
# combined_data <- SCTransform(combined_data, verbose = TRUE,batch_var="ltx_case",conserve.memory=TRUE,method = 'glmGamPoi',residual_type="pearson",min_cells=10)
# combined_data <- SCTransform(combined_data, verbose = TRUE,batch_var="ltx_case",conserve.memory=TRUE,min_cells=10)
combined_data <- SCTransform(combined_data, verbose = TRUE,conserve.memory=TRUE,min_cells=10)


combined_data<- FindVariableFeatures(combined_data, selection.method = "vst", nfeatures = 2000)
#perform linear reduction
combined_data <- RunPCA(combined_data, features = VariableFeatures(object = combined_data))
pdf(gdpath("pc_loadings.pdf",root=subroot,external=root_external))
print(VizDimLoadings(combined_data, dims = 1:2, reduction = "pca"))
print(DimPlot(combined_data, reduction = "pca"))
DimHeatmap(combined_data, dims = 1:15, cells = 500, balanced = TRUE)
#JackStraw Disables by Seurat authors as it has been found it is not good for SCTransformed data
# combined_data <- JackStraw(combined_data, num.replicate = 100)
# combined_data <- ScoreJackStraw(combined_data, dims = 1:pc_use)
# print(JackStrawPlot(combined_data, dims = 1:20))
print(ElbowPlot(combined_data,ndims=30))
dev.off()

combined_data <- FindNeighbors(combined_data, dims = 1:pc_use)
combined_data <- FindClusters(combined_data, resolution = cluster_resolution)

pdf(gdpath("pC_norm.pdf",root=subroot,external=root_external))
combined_data<-RunUMAP(combined_data,dims=1:pc_use)
print(DimPlot(combined_data, reduction = "umap",group.by=c("timepoint","seurat_clusters")))
dev.off()
cat("Simple per-cellnormalization complete.\n")
mcsaveRDS(combined_data,file=gdpath(paste0("combined_data_no_anchored_data",".rds"),root=subroot,external=root_external))
rm(combined_data)
gc()
###################################################################################################
# End of the simple per-cell normalization technique
###################################################################################################


