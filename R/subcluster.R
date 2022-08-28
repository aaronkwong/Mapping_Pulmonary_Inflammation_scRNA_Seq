# this script reads a tab delimted file with where each row specifies a cluster for subclustering. There are two columns, the first column specifies the name of the cluster which should be subclustered, and the second column specifies the number of sub communities which should be found in that cluster. This script outputs a new seurat_object with the subclustering applied. 
library(renv)
renv::activate()

library(trqwe)
library(Seurat)
library(patchwork)

source("gdpath.R")
source("subcluster_functions.R")


suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-r", "--pc_use"), default="NA",
            help="The number of principal components to use when subclustering to find new sub communities"),
  #
  make_option(c("-c", "--project_name"), default="NA",
            help="A name for the project."),
  #
  make_option(c("-a", "--input_dir"), default="NA",
            help="The directory containing the combined and anchored rds (within the root directory)"),  
  #
  make_option(c("-v", "--root"), default="NA",
        	help="The directory where the various result directories for each part of the pipeline will be written to."),
  #
  make_option(c("-q", "--root_external"), default="NA",
            help="TRUE or FALSE. If the root is contatined in the cloud synced folder."),
  #
  make_option(c("-b", "--make_save"), default=FALSE,
            help="TRUE or FALASE. Should subclustered seurat object be saved?")
)

opt <- parse_args(OptionParser(option_list=option_list))



#read in command line variables
n_pcs<-as.numeric(opt$pc_use)
project_name<-opt$project_name
root<-opt$root
input_dir<-opt$input_dir

if(opt$make_save=="TRUE"){
	make_save<-TRUE
}else{
	make_save<-FALSE
}

if(opt$root_external=="TRUE"){
  root_external<-TRUE
}else{
  root_external<-FALSE
  cat("internal root detected.\n")
}

#paths to saved files for processing
combined_path<-paste0(root,input_dir,project_name,"_combined_data.rds")
anchored_path<-paste0(root,input_dir,project_name,"_ANCHORED_integration.rds")

#create root for output directory
subroot<-paste0(root,"subclustering/")
dir.create(subroot)

#check to see if input file specifies any clusters that should be sub clustered
exit_value<-0
while(exit_value!=1){
	#check if the subclustering file which specifies clusters for subclustering exists
	if(file.exists(paste0(root,"subclustering/subclustering_data.csv"))){
		cat("found subclutering data..\n")
		subclustering_data<-read.csv(paste0(root,"subclustering/subclustering_data.csv"))
		print(subclustering_data)
		#if file exists but no information in it, then just copy input to output
		if(nrow(subclustering_data)==0){
			cat("no data found in subclustering from...\n")
			anchored_data<-mcreadRDS(anchored_path)
			mcsaveRDS(anchored_data,gdpath("project_data_ANCHORED_integration.rds",root=subroot,external=root_external),mc.cores=16)
			combined_data<-mcreadRDS(combined_path)
			mcsaveRDS(combined_data,gdpath("combined_data.rds",root=subroot,external=root_external),mc.cores=16)
			quit()
		}else{
			exit_value<-1
		}
	}else{
		#if no file specifying clusters for subclustering then make one
		cluster_for_subclustering<-c(as.numeric())
		expected_communities<-c(as.numeric())
		subclustering_data<-data.frame(cluster_for_subclustering,expected_communities)
		write.table(subclustering_data,gdpath("subclustering_data.csv",root=subroot,external=root_external),row.names=FALSE,sep=",")
	}
}

#read in data
cat("reading in anchored data...\n")
anchored_data<-mcreadRDS(anchored_path)

# plot clusters before doing subclustering of cells in anchored space
pdf(gdpath("subclustering_plots.pdf",root=subroot,external=root_external))
	DimPlot(anchored_data,label=TRUE)+ plot_annotation(title = paste0("anchored_original"))


	# for each row specific in the subclustering table, subcluster to the desired number of sub populations
	for (i in 1:nrow(subclustering_data)){
		anchored_data<-subcluster(anchored_data=anchored_data,cluster_for_subclustering=subclustering_data$cluster_for_subclustering[i],desired_subclusters=subclustering_data$expected_communities[i],n_pcs=n_pcs)
	}

	#save anchored object with subclustered labels
	Idents(anchored_data)<-anchored_data$integrated_clusters
	cat("saving modified anchored data...\n")
	if(make_save){
		mcsaveRDS(anchored_data,gdpath("project_data_ANCHORED_integration.rds",root=subroot,external=root_external),mc.cores=16)
	}

	# plot clusters after subclustering of cells in anchored space
	DimPlot(anchored_data,label=TRUE)+ plot_annotation(title = paste0("anchored_post_subclustering"))

	# read in data and apply the subclustered labels
	cat("reading in combined data...\n")
	combined_data<-mcreadRDS(combined_path)
	# plot clusters before doing subclustering of cells in combined space
	DimPlot(combined_data,label=TRUE)+ plot_annotation(title = paste0("combined_original"))
	combined_data<-AddMetaData(combined_data,metadata = anchored_data$integrated_clusters,col.name="integrated_clusters")
	Idents(combined_data)<-combined_data$integrated_clusters
	cat("saving modified combined data...\n")
	if(make_save){
		mcsaveRDS(combined_data,gdpath("combined_data.rds",root=subroot,external=root_external),mc.cores=16)
	}

	# plot clusters after subclustering of cells in anchored space
	DimPlot(combined_data,label=TRUE)+ plot_annotation(title = paste0("combined_post_subclustering"))
dev.off()