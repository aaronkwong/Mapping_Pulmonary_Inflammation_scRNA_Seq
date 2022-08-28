suppressPackageStartupMessages(library(optparse))
library(Seurat)
library(trqwe)
source("gdpath.R")
suppressPackageStartupMessages(library(future))
options(future.globals.maxSize = 1200000 * 1024^2)
plan("multiprocess", workers = 16)

option_list <- list(
  # make_option(c("-h", "--help"), action="store_true", default=FALSE, 
  #              help="Show this help message and exit"),
  
  make_option(c("-r", "--root"), default="NA",
              help=""),
  #
  make_option(c("-t", "--root_external"), default="NA",
              help=""),  
  #
  make_option(c("-a", "--pc_use"), default="NA",
              help=""),
  #
  make_option(c("-b", "--project_name"), default="NA",
              help=""),
  #
  make_option(c("-d", "--project_data_path"), default="NA",
              help=""),
  #
  make_option(c("-c", "--cluster_resolution"), default="NA",
              help=""),
  #
  make_option(c("-e", "--remove_genes"), default="NA",
              help=""),
  #
  make_option(c("-f", "--remove_doublets"), default="TRUE",
              help="")
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
remove_genes<-unlist(strsplit(opt$remove_genes,","))

source(gdpath("Masters/project_6_scRNA_seq/scRNA_seq_Cell_Tracker/combined_analysis/integration_functions.R"))
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

# ###################################################################################################
# # Start of the Anchored Integrated  normalization technique rPCA
# ###################################################################################################
# cat("Starting Anchored Integration of two timepoints....","\n")
# #combine two datasets into a single object
# # combined_data_list <- c(pre.evlp,post.evlp)

# #normalize each dataset
# combined_data_list<-list()
# for (i in 1:length(all_samples)) {
#   combined_data_list[[i]] <- SCTransform(all_samples[[i]], verbose = T)
# }

# #remove all_samples to free up ram as we no longer need it
# rm(all_samples)
# gc()

# #rpca may be needed here due to memory limitations https://satijalab.org/seurat/articles/integration_large_datasets.html
# #calculate pearson features
# combined_data.features <- SelectIntegrationFeatures(object.list = combined_data_list, nfeatures = 3000)

# #un PCA for rpca
# combined_data_list<-lapply(X = combined_data_list, FUN = function(x) {
#     # x <- ScaleData(x, features = features, verbose = FALSE) #can remove RunPCA pulls from scale.data slot of SCT
#     x <- RunPCA(x, features = combined_data.features, verbose = FALSE)
# })


# combined_data_list <- PrepSCTIntegration(object.list = combined_data_list, anchor.features = combined_data.features, verbose = FALSE)

# #calculate anchors
# # combined_data.anchors <- FindIntegrationAnchors(object.list = combined_data_list, normalization.method = "SCT", 
# #   anchor.features = combined_data.features, verbose = TRUE)

# combined_data.anchors <- FindIntegrationAnchors(object.list = combined_data_list,anchor.features=combined_data.features, normalization.method = "SCT", verbose = TRUE,reduction="rpca")
# combined_data.integrated <- IntegrateData(anchorset = combined_data.anchors, normalization.method = "SCT", 
#   verbose = TRUE)

# #display umap
# combined_data.integrated <- RunPCA(combined_data.integrated, verbose = FALSE)
# combined_data.integrated <- RunUMAP(combined_data.integrated, dims = 1:30)
# combined_data.integrated <- FindNeighbors(combined_data.integrated, dims = 1:pc_use)
# combined_data.integrated <- FindClusters(combined_data.integrated, resolution = cluster_resolution)
# #this saves the anchored object with the default assay as anchored. This should only be used cluster label annotation.
# mcsaveRDS(combined_data.integrated,file=gdpath(paste0("project_data","_ANCHORED_integration",".rds"),root=subroot,external=root_external))

# cat("Anchored Integration of two timepoints complete.","\n")
# ###################################################################################################
# # End of the Anchored Integrated  normalization technique
# ###################################################################################################

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
plots <- DimPlot(combined_data.integrated, group.by = c("seurat_clusters","timepoint"), combine = FALSE,label=T)
print(CombinePlots(plots))
plots <- DimPlot(combined_data.integrated, group.by = c("ltx_case"), combine = FALSE,label=T)
print(CombinePlots(plots))
plots <- DimPlot(combined_data.integrated, group.by = c("seurat_clusters"), combine = FALSE,label=T)
print(CombinePlots(plots))
dev.off()
