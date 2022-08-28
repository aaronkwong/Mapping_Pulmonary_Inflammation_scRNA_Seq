# This script loads single cell data from each case and combines data into a single object. The script outputs a single combined seurat object and some summary plots.

suppressPackageStartupMessages(library(renv))
renv::activate()
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sctransform))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(scClustViz))
suppressPackageStartupMessages(library(trqwe))
suppressPackageStartupMessages(library(future))

#set analysis options
options(future.globals.maxSize = 120000 * 1024^2)
plan("multiprocess", workers = 16)

#source custom functions
source("gdpath.R")
source("integration_functions.R")

option_list <- list(
    make_option(c("-r", "--root"), default="NA",
              help="The directory where the various result directories for each part of the pipeline will be written too."),
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
    make_option(c("-e", "--project_directory"), default="NA",
              help="Directory containing the quality controlled files."),
    #
    make_option(c("-x", "--max_cell_scaling"), default=0,
              help="The maximum number of cells for analysis. If total cells exceeds this number, each sample will be randomly scaled down to be under this max.")
    )

opt <- parse_args(OptionParser(option_list=option_list))


#seed for reproduceability
set.seed(42)

#set the project directory
project_directory<-opt$project_directory
#set the directory containing project data file
project_data_path<-opt$project_data_path
root<-opt$root
#convert text input to logical
if(opt$root_external=="TRUE"){
    root_external<-TRUE
}else if(opt$root_external=="FALSE"){
    root_external<-FALSE
}else{
  cat("root_external was not defined as true or flase.\n")
}

#set clustering resolution
cluster_resolution=as.numeric(opt$cluster_resolution)
#number of PC's to use 
pc_use=as.numeric(opt$pc_use)
#set the maximum number of cells to allow in the analysis. Will downsample the number of cells from each sample proportionately if this number is exceeded
max_cell_scaling=as.numeric(opt$max_cell_scaling)

#read in the project file
cat(paste0("reading in project_data at: ",project_data_path,"\n"))
project_data<-read.csv(project_data_path)
# only take project data rows with a directory because some may be empty
project_working<-project_data[project_data$directory!="",]
cat(paste0("Found ",nrow(project_data)," potential samples found in the project data. Of those ",length(unique(project_working$case))," are populated.\n"))
cat(paste0("The following cases will be used in the analysis: \n",paste0(project_working$directory,collapse="\n"),"\n"))

#set the output directory in the results folder
subroot<-paste0(root,"initial_integration/")
dir.create(subroot)

#lets read in the quality controlled data for each case
#first lets get a list of all the files we need to load
files_to_load<-c()
for (i in 1:nrow(project_working)){
	#put together the path for the results directory for samples
	files_to_load<-c(files_to_load,paste0(project_directory,project_working$directory[i],"/quality_control/",project_working$timepoint[i],"/",project_working$sample_name[i],"_qc.rds"))
}
project_working$qc_path<-files_to_load

#now lets load files into memory
all_samples<-list()
for (i in 1:nrow(project_working)){
	all_samples[[i]]<-mcreadRDS(project_working$qc_path[i],mc.cores=16)
  # read in data
  cat(paste0("Read in data: ",project_working$qc_path[i],". Found ",ncol(all_samples[[i]]), " cells in this sample.\n"))

	# assign a meta data defining whether pre or post
	if (project_working$timepoint[i]=="pre"){
		cluster_letters <- rep("pre",length(Idents(object = all_samples[[i]])))
	}else if(project_working$timepoint[i]=="post"){
		cluster_letters <- rep("post",length(Idents(object = all_samples[[i]])))
	}
    names(cluster_letters) <- colnames(x = all_samples[[i]])
    all_samples[[i]] <- AddMetaData(
	    object = all_samples[[i]],
	    metadata = cluster_letters,
	    col.name = 'timepoint'
    )

    #add the case from which the sample came
    ltx_case_label <- rep(project_working$case[i],length(Idents(object = all_samples[[i]])))
    names(ltx_case_label) <- colnames(x = all_samples[[i]])
    all_samples[[i]] <- AddMetaData(
	    object = all_samples[[i]],
	    metadata = ltx_case_label,
	    col.name = 'ltx_case'
    )

    #add the sample_name from which the sample came
    sample_name_label <- rep(project_working$sample_name[i],length(Idents(object = all_samples[[i]])))
    names(sample_name_label) <- colnames(x = all_samples[[i]])
    all_samples[[i]] <- AddMetaData(
	    object = all_samples[[i]],
	    metadata = sample_name_label,
	    col.name = 'sample_name'
    )

    #because different samples can have cells that use the same barcode we need to append an ending to each sample so that duplicate cell names dont collide
    all_samples[[i]] <- RenameCells(object = all_samples[[i]], add.cell.id = paste0(project_working$case[i],project_working$timepoint[i]))
}

# if max_cell_scaling does not equal zero (0 is defualt and will use all cells)
if(max_cell_scaling!=0){
    cat("\nmax_cell_scaling selected. Scaling down samples:\n")
    total_cells<-0
    for (i in 1:length(all_samples)){
        total_cells<-total_cells+ncol(all_samples[[i]])
    }
    cat(paste0("there are a total of ",total_cells," in this analysis.\n"))
    if(total_cells > max_cell_scaling){
        for (i in 1:length(all_samples)){
            num_cells_in_sample<-ncol(all_samples[[i]])
            index_cells_to_use<-sample(x=1:num_cells_in_sample,size=round((num_cells_in_sample/total_cells)*max_cell_scaling),replace=FALSE)
            all_samples[[i]]<-(all_samples[[i]])[,index_cells_to_use]
            cat(paste0(num_cells_in_sample," --> ",length(index_cells_to_use),"\n"))
        }
    }else{
        cat(paste0("max cell scaling selected, but the total number of cells is less than ",max_cell_scaling," cells as set by max_cell_scaling\n"))
    }
}


###################################################################################################
# Start of the simple per-cell normalization technique
###################################################################################################
cat("Starting simple per-cellnormalization...\n")
#combine two datasets into a single object
combined_data <- merge(x=all_samples[[1]],y=all_samples[-1])
#set min_cells=10 or vst can throw errors
combined_data <- SCTransform(combined_data, verbose = TRUE,conserve.memory=TRUE,min_cells=10)
#find variable features
combined_data<- FindVariableFeatures(combined_data, selection.method = "vst", nfeatures = 2000)
#perform linear reduction
combined_data <- RunPCA(combined_data, features = VariableFeatures(object = combined_data))
#print out QC metrics
pdf(gdpath("pc_loadings.pdf",root=subroot,external=root_external))
    print(VizDimLoadings(combined_data, dims = 1:2, reduction = "pca"))
    print(DimPlot(combined_data, reduction = "pca"))
    DimHeatmap(combined_data, dims = 1:15, cells = 500, balanced = TRUE)
    print(ElbowPlot(combined_data,ndims=30))
dev.off()

#find neighbours for clustering
combined_data <- FindNeighbors(combined_data, dims = 1:pc_use)
combined_data <- FindClusters(combined_data, resolution = cluster_resolution)

#construct preliminary UMAP
pdf(gdpath("pC_norm.pdf",root=subroot,external=root_external))
    combined_data<-RunUMAP(combined_data,dims=1:pc_use)
    print(DimPlot(combined_data, reduction = "umap",group.by=c("timepoint","seurat_clusters")))
dev.off()

cat("Simple per-cellnormalization complete.\n")
#save results
mcsaveRDS(combined_data,file=gdpath(paste0("combined_data_no_anchored_data",".rds"),root=subroot,external=root_external),mc.cores=16)
###################################################################################################
# End of the simple per-cell normalization technique
###################################################################################################
