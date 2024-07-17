suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-t", "--seurat_object_path_with_recip_cells"), default="NA",
        help="a seurat object with the recipient cells. This algo will seperate out the recip cells from meta data and give them a seperate cluster"),
  #
    make_option(c("-x", "--cluster_resolution"), default="NA",
        help=""),
  #
    make_option(c("-y", "--include_recip"), default="NA",
        help=""),
	make_option(c("-b", "--clusters_containing_recip_cells"), default="NA",
        help="MAKE SURE clusters are in ascending order and this matches the manual anno with recip in the manual_changes folder"),
    make_option(c("-a", "--root"), default="NA",
    	help=""),   
  #
    make_option(c("-p", "--root_external"), default="NA",
        help=""),
  #
    make_option(c("-c", "--project_name"), default="NA",
        help="")
  )

opt <- parse_args(OptionParser(option_list=option_list))

library(CellChat)
library(trqwe)
library(Seurat)
future::plan("multiprocess", workers = 16) # do parallel
# options(future.fork.multithreading.enable = FALSE)
source("gdpath.R")

#read vars
root<-opt$root
cluster_resolution<-opt$cluster_resolution
project_name<-opt$project_name
seurat_object_path_with_recip_cells<-opt$seurat_object_path_with_recip_cells
clusters_expected_to_contain_recip_cells<-as.numeric(unlist(strsplit(x=opt$clusters_containing_recip_cells,split=",")))
if(any(is.na(clusters_expected_to_contain_recip_cells))){
	cat("error. one of the clusters input as containing recipient cells is not an integer.\n")
	quit()
}

#check whether recipient cells should be included in analysis
if (opt$include_recip=="TRUE"){
  include_recip<-TRUE
}else{
  include_recip<-FALSE
}

#convert text input to logical for root_external
if(opt$root_external=="TRUE"){
    root_external<-TRUE
}else if(opt$root_external=="FALSE"){
    root_external<-FALSE
}else{
  cat("root_external was not defined as true or flase.\n")
}

#make result subroot
subroot<-paste0(root,"cellchat/")
dir.create(subroot)

#read in the seurat object, which should be the QC (doublet removed object)
seurat_object<-mcreadRDS(seurat_object_path_with_recip_cells)


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

seurat_object_with_recipient<-seurat_object[,seurat_object$timepoint=="post"]
seurat_object_pre<-seurat_object[,seurat_object$timepoint=="pre"]

# we need to clean cells. There may be some clusters that contain a single called recipient cell. We need to specify the cluster in which
# we are expecting recipient cells. Clusters that we dont have expect recipient cells in should have those "recipient" cells removed.
cat(paste0("Removing recipient cells from clusters they are not expected to be in. Keeping recip cells only in cluster: ",clusters_expected_to_contain_recip_cells,"\n"))
seurat_object_with_recipient_clean<-seurat_object_with_recipient[,!((!(seurat_object_with_recipient@meta.data$integrated_clusters %in% clusters_expected_to_contain_recip_cells)) & (seurat_object_with_recipient@meta.data$recipient_origin==TRUE))]
table(seurat_object_with_recipient_clean@meta.data$integrated_clusters)
#next we should give each cluster of recipient cells a unique cluster in the integrated clusters column of the metadata
#get the biggest cluster number FROM THE WHOLE DATASET. WHAAT IF THE MAX CLUSTER EXISTS ONLY IN PRE TIMEPOINT
new_cluster_label<-max(as.numeric(na.omit(levels(seurat_object_with_recipient@meta.data$integrated_clusters))))
#create a new vector of labels and modify recipient cells as their own seperate cluster
integrated_clusters<-as.numeric(levels(seurat_object_with_recipient_clean@meta.data$integrated_clusters)[seurat_object_with_recipient_clean@meta.data$integrated_clusters])
sink(gdpath("new_recip_cluster.txt",root=subroot,external=root_external))
for (x in clusters_expected_to_contain_recip_cells){
	#add 1 to the current largest cluster as the new label we will use. This will increase by 1 every loop if there is more than 1 cluster containing recipient cells
	new_cluster_label<-new_cluster_label+1
	cat(paste0("deconvolving cluster: ",x,". New cluster label: ",new_cluster_label,"\n"))
	#replace recipient cells label with new label number
	print(sum(seurat_object_with_recipient_clean@meta.data$integrated_clusters==x & seurat_object_with_recipient_clean@meta.data$recipient_origin==TRUE))
	if(sum(seurat_object_with_recipient_clean@meta.data$integrated_clusters==x & seurat_object_with_recipient_clean@meta.data$recipient_origin==TRUE)==0){
		cat("Error, one of the clusters supposed to contain recipient cells did not have any.\n")
		quit()
	}
	integrated_clusters[seurat_object_with_recipient_clean@meta.data$integrated_clusters==x & seurat_object_with_recipient_clean@meta.data$recipient_origin==TRUE]<-new_cluster_label
}
sink()

#replace integrated cluster names with the manual names with recipient labels file if it exists
manual_anno_with_dr<-paste0(root,"manual_changes/manual_cluster_names_with_recip_clusters.tab")
if(file.exists(manual_anno_with_dr)){
	man_anno<-read.csv(manual_anno_with_dr,stringsAsFactors=FALSE,sep="\t")
	if(length(man_anno$integrated_clusters)==length(man_anno$manual)){
		cat(paste0("cluster annotations are acceptable."))
		for(i in 1:nrow(man_anno)){
			integrated_clusters[integrated_clusters==man_anno$integrated_clusters[i]]<-man_anno$manual[i]
		}
	}
}


#factor the new labels
#use new factored labels and overwrite old integrated labels
seurat_object_with_recipient_clean<-AddMetaData(object=seurat_object_with_recipient_clean, metadata=integrated_clusters, col.name = "integrated_clusters")
seurat_object_with_recipient_clean@meta.data$integrated_clusters<-integrated_clusters


# we need to redo the labels on the 
integrated_clusters_pre<-as.numeric(levels(seurat_object_pre@meta.data$integrated_clusters)[seurat_object_pre@meta.data$integrated_clusters])
manual_anno_with_dr<-paste0(root,"manual_changes/manual_cluster_names_with_recip_clusters.tab")
if(file.exists(manual_anno_with_dr)){
	man_anno<-read.csv(manual_anno_with_dr,stringsAsFactors=FALSE,sep="\t")
	if(length(man_anno$integrated_clusters)==length(man_anno$manual)){
		cat(paste0("cluster annotations are acceptable."))
		for(i in 1:nrow(man_anno)){
			integrated_clusters_pre[integrated_clusters_pre==man_anno$integrated_clusters[i]]<-man_anno$manual[i]
		}
	}
}
seurat_object_pre<-AddMetaData(object=seurat_object_pre, metadata=integrated_clusters_pre, col.name = "integrated_clusters")

run_cell_chat<-function(seurat_object){
	#log normalized data with pseudocount of 1 (this is log1p)
	data.input <- seurat_object@assays$SCT@data
	#extract cell cluster labels
	meta = data.frame(integrated_clusters=paste0(seurat_object$integrated_clusters))
	rownames(meta)<-names(seurat_object$integrated_clusters)

	#list all the labels
	unique(meta$integrated_clusters)

	#create the cellchat object
	cellchat <- createCellChat(object = data.input, meta = meta, group.by = "integrated_clusters")

	#set the cell chat database to use
	CellChatDB <- CellChatDB.human
	CellChatDB.use <- CellChatDB # simply use the default CellChatDB

	# set the used database in the object. Use the full database
	cellchat@DB <- CellChatDB.use

	# subset the expression data of signaling genes for saving computation cost
	cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

	cellchat <- identifyOverExpressedGenes(cellchat)
	cellchat <- identifyOverExpressedInteractions(cellchat)
	# project gene expression data onto PPI network (optional)
	cellchat <- projectData(cellchat, PPI.human)

	cellchat <- computeCommunProb(cellchat)
	# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
	cellchat <- filterCommunication(cellchat, min.cells = 10)

	#get pathway level interactions
	cellchat <- computeCommunProbPathway(cellchat)

	#LR interaction network
	cellchat <- aggregateNet(cellchat)

	return(cellchat)
}

mcsaveRDS(run_cell_chat(seurat_object=seurat_object_with_recipient_clean),file=gdpath(paste0(project_name,"_post_with_recipient_cellchat.rds"),root=subroot,external=root_external))
mcsaveRDS(run_cell_chat(seurat_object=seurat_object_pre),file=gdpath(paste0(project_name,"_pre_cellchat.rds"),root=subroot,external=root_external))
