# This script takes an integrated seurat object (and the clustering information) and identifies doublets as described in "DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors" by Christopher McGinnis. Outputs a seurat object with an additional metadata columns specifying the doublet score and doublet prediction. 

library(renv)
renv::activate()
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(trqwe))
suppressPackageStartupMessages(library(DoubletFinder))

#source custom functions
source("gdpath.R")
source("doubletFinder_functions.R")


option_list <- list(
  make_option(c("-r", "--root"), default="NA",
              help="The directory where the various result directories for each part of the pipeline will be written to."),
  #
  make_option(c("-t", "--root_external"), default="NA",
              help="TRUE or FALSE. If the root is contatined in the cloud synced folder."),  
  #
  make_option(c("-a", "--pc_use"), default="NA",
              help="Number of PCs to use when clustering artifical droplets. "),
  #
  make_option(c("-b", "--project_name"), default="NA",
              help="A name for the project."),
  #
  make_option(c("-c", "--project_data_path"), default="NA",
              help="Directory containing the quality controlled files.")
  )

opt <- parse_args(OptionParser(option_list=option_list))

options(future.globals.maxSize = 120000 * 1024^2)
suppressPackageStartupMessages(library(future))
plan("multiprocess", workers = 16)


#load in arguments
root<-opt$root
#convert text input to logical
if(opt$root_external=="TRUE"){
    root_external<-TRUE
}else if(opt$root_external=="FALSE"){
    root_external<-FALSE
}
#project name
project_name<-opt$project_name
#pcs to use
pc_use<-as.numeric(opt$pc_use)
#read in project data and valid rows
project_data_path<-opt$project_data_path
project_data<-read.csv(project_data_path)
project_working<-project_data[project_data$directory!="",]


#read in saved data file 
cat("reading in data files.\n")
combined_data<-mcreadRDS(gdpath(paste0("initial_integration/combined_data",".rds"),root=root,external=root_external))

#create the subdirectory for results
subroot<-paste0(root,"doublet_removal/")
dir.create(subroot)


# This script will make doublet finder scans through a range of QC parameters on the first run
# and after analysis lets the user choose which parameters they think works best. Once prefered 
# parameters are selected they should be written into the "DF_parameters.csv" file. Here the script
# Checks to see this file exists, and if so to read the parameters in for doublet removal.
if (file.exists(paste0(root,"doublet_removal/DF_parameters.csv"))){
	cat("previous DF_parameter file found, reading it in.\n")
	df_parameters<-read.csv(paste0(root,"doublet_removal/DF_parameters.csv"),row.names=1,header=TRUE)
	print(df_parameters)
	if (!all(project_working$sample_name==rownames(df_parameters))){
		cat("The sample names in project_working$sample_name did not match rownames(df_parameters)\n")
		quit()
	}
# if no "DF_parameters.csv" exists, this is the first run. Scan a range of parameters for every 
# single cell sample.
# save the results for review
}else{ 
	first_run<-TRUE
	cat("first run detected. Creating DF_parameter file.\n")
	pK<-rep(0,length(project_working$sample_name))
	pN<-rep(0,length(project_working$sample_name))
	poisson_multiplet_rate<-rep(0,length(project_working$sample_name))
	df_parameters<-data.frame(pK,pN,poisson_multiplet_rate)
	rownames(df_parameters)<-project_working$sample_name
	write.table(df_parameters,file=gdpath("DF_parameters.csv",root=subroot,external=root_external),sep=",",col.name=NA)
}

#list the batches
cat("list of batches that will have doublets removed:\n")
cat(paste0("subsetting cells from ",project_working$case," ",project_working$timepoint,".\n"))

# check to see if user selected doublet removal numbers have been entered for all samples
# If at least one row is all zeros then we need to run a sweep 
# to find the best parameters for that sample and it is added to the que "samples_to_process"
# 
# If no rows are all 0, then this means that all rows must usser specific paramaters entered and we can proceed with 
# the actual doublet removal. In this case all samples are added to the samples_to_process que
#
# The que is smart enough to know how to whether an item in que needs the sweep or is ready 
# for doublet removal with user specific parameters
if(any(apply(df_parameters,1,FUN=function(x){all(x==0)}))){ # if at least one row of all zeros found
	samples_to_process<-which(apply(df_parameters,1,FUN=function(x){all(x==0)}))
	cat(paste0("Found samples requiring sweep. Sweep will be conducted on samples: ",paste0(rownames(df_parameters)[samples_to_process],collapse=","),"\n"))
	#because the df parameters file exists it wont be considered a first run, even though it really is because samples require sweep. Set this value back to true
	first_run<-TRUE
}else{
	samples_to_process<-1:length(unique(rownames(df_parameters)))
	cat(paste0("all samples have a value entered for doublet removal. Proceesing with doublet removal.\n"))
	first_run<-FALSE
}

# for each sample in the que, pull the appropriate cells from the seurat object and process
cat("running doublet analysis ...\n")
combined_data_doublet_finder<-list()
x<-1
for (i in samples_to_process){
	#subset big object into each case by pre and post
	cat(paste0("subsetting cells from ",project_working$case[i]," ",project_working$timepoint[i],".\n"))
	temp_subset<-combined_data[,combined_data@meta.data$ltx_case==project_working$case[i] & combined_data@meta.data$timepoint==project_working$timepoint[i]]
	# if the pN and pK values passed in are 0 (user has yet to specify) then a parameter sweep is 
	# conducted. If non zero values, doublet removal is conducted. 
	combined_data_doublet_finder[[x]]<-run_find_doublets(
	    seurat_object_processed=temp_subset,
	    pcs=pc_use,
	    poission_doublet_rate=df_parameters$poisson_multiplet_rate[i],
	    project_name=project_name,
	    sample_name=rownames(df_parameters)[i],
	    manually_selected_pN=df_parameters$pN[i],
	    manually_selected_pK=df_parameters$pK[i],
	    root=subroot,
	    root_external=TRUE
	)
	x<-x+1
	cat(paste0("\nadding results. ",length(combined_data_doublet_finder),".\n"))
}

#check to see if the output of fun_find_doublets is a seurat object. If this is the first run that ran sweeps there would be nothing to save
if (!first_run){
	#now lets combine all seperated seurat objects
	combined_data_doublet_finder<-merge(combined_data_doublet_finder[[1]],combined_data_doublet_finder[-1])

	#every time doublet finder runs on each sample it adds 2 new meta data columns for doublet score and classification
	#this means if we have 4 scRNA seq samples in our dataset we will have 8 extra columns generated
	#ultimately we want just a single metadata column with the score and doublet classification
	#we need combine all doublet scores into single meta data column. Also do the same for classification
	#find all columns with score
	doublet_scores<-combined_data_doublet_finder@meta.data[,grepl("pANN",colnames(combined_data_doublet_finder@meta.data))]
	#find all columns with classficiation
	doublet_classification<-combined_data_doublet_finder@meta.data[,grepl("DF.classifications",colnames(combined_data_doublet_finder@meta.data))]
	#collapse all scores into 1 column
	doublet_scores_clean<-apply(doublet_scores,MARGIN=1,FUN=function(x){
		if(length(na.omit(x)==1)){
			return(na.omit(x))
		}else{
			cat("error, one of the cells has two doublet scores. was a cell scored twice?")
			quit()
		}
	})
	#collapse all classifcation into one column
	doublet_classification_filtered<-apply(doublet_classification,MARGIN=1,FUN=function(x){
		if(length(na.omit(x)==1)){
			return(na.omit(x))
		}else{
			cat("error, one of the cells has two doublet scores. was a cell scored twice?")
			quit()
		}
	})
	# add scores and doublet classification meta data 
	combined_data_doublet_finder <- AddMetaData(object = combined_data_doublet_finder, metadata = doublet_scores_clean, col.name = 'DoubletFinderScores')
	combined_data_doublet_finder <- AddMetaData(object = combined_data_doublet_finder, metadata = doublet_classification_filtered, col.name = 'DoubletFinderClassification')

	#make a quick plot to refer back to
	pdf(gdpath("hist_doublet_scores_by_case.pdf",root=subroot,external=root_external))
	for (case in unique(combined_data_doublet_finder$ltx_case)){
		hist(combined_data_doublet_finder$DoubletFinderScores[combined_data_doublet_finder$ltx_case==case & combined_data_doublet_finder$timepoint=="pre"])
		hist(combined_data_doublet_finder$DoubletFinderScores[combined_data_doublet_finder$ltx_case==case & combined_data_doublet_finder$timepoint=="pre"])
	}
	dev.off()


	#save doublet result for future review
	mcsaveRDS(combined_data_doublet_finder,file=gdpath("combined_data_doublet_finder.rds",root=subroot,external=root_external))
}

#sink a summary of the new Seurat object
sink(gdpath("combined_data_summary.txt",root=subroot,external=root_external))
combined_data_doublet_finder
print(table(combined_data_doublet_finder$DoubletFinderClassification,combined_data_doublet_finder$ltx_case))
sink()


