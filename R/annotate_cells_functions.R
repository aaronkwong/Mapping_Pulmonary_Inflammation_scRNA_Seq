library(Seurat)
library(stringr)

#this function extracts the cells from the pre timepoint
extract_cells_using_timepoint<-function(seurat_object,timepoint){
	pre_cells<-seurat_object[,seurat_object$timepoint==timepoint]
	return(pre_cells)
}

create_marker_table<-function(seurat_object,timepoint="pre",min_prct_cells,logfc.threshold,root,root_external){
	#get the cells from only the timpoint we want, usually this is pre
	seurat_object<-extract_cells_using_timepoint(seurat_object,timepoint)
	marker_out <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = min_prct_cells, logfc.threshold = logfc.threshold)
	#we could further cutoff the marker_put file by an fdr, but even if no markers meet cutoff for a given cell type I think
	#we would still want just the top markers anyways regardless of fdr. So lets not cut off
	#get a list of the unique clusters
	unique_clusters<-levels(Idents(seurat_object))
	marker_table<-list()
	i<-1
	for (clus in unique_clusters){
		marker_table[[i]]<-marker_out[marker_out[,"cluster"]==clus,"gene"]
		i<-i+1
	}
	#this functions pads the list so a table can be written
	pad_list_col<-function(x){
	 	max_row_size<-max(sapply(x,length))
	 	new_table<-list()
	 	for (i in 1:length(x)){
	 		new_table[[i]]<-c(x[[i]],rep(NA,times=max_row_size-length(x[[i]])))
	 	}
	 	new_table<-do.call(cbind,new_table)
	 	new_table<data.frame(new_table)
	 	rownames(new_table)<-names(x)
	 	return(new_table)
	 }
	fin<-pad_list_col(marker_table)
	colnames(fin)<-unique_clusters
	#write out result
	write.table(fin,gdpath(paste0("marker_table_",timepoint,"_only_nonint.tab"),root=root,external=root_external),sep="\t",row.names=F)
}

annotate_combined_using_pre_only<-function(seurat_object,infile_gmt,root,root_external){
	#subset the seurat object into cells from the pre timepoint only
	seurat_object_pre<-extract_cells_using_timepoint(seurat_object,timepoint="pre")
	seurat_object_post<-extract_cells_using_timepoint(seurat_object,timepoint="post")
	#we need to make sure that all the clusters are found in both timepoints. If there is a cluster only present in the post timepoint
	#then we need to use a combination of both
	#if idents from pre and post are the same we will use only the pre timepoint for automatic annotation
	if (identical(seurat_object_pre,seurat_object_post)){
		seurat_object<-seurat_object_pre
		cat("Pre timepoint used.\n")
	#if there are clusters that dont exist in both timepoints then we will just use all cells from the pre and post timepoints
	}else{
		sink(gdpath("both_timepoints_used_for_integration.txt",root=root,external=root_external))
		cat("TRUE")
		sink()
		if(ncol(seurat_object)>80000){
			cat("too many cells, subsampling to 80,000 cells for GSVA annotation only\n")
			seurat_object<-seurat_object[,sample(1:ncol(seurat_object),80000)]
		}else{
			seurat_object<-seurat_object
		}
		cat("Combined Pre and Post timepoints used.\n")
	}
	#AverageExpression will trip an error if there are any NA values in expression. Genes that dont meet min.cells in SCTransform will all be NA. must remove
	seurat_object_average_expression<-seurat_object
	seurat_object_average_expression@assays$SCT@data<-seurat_object@assays$SCT@data[!apply(seurat_object@assays$SCT@data,1,FUN=function(x){all(is.na(x))}),]
	#use cleaned seurat object on AverageExpression
    avg_expression<-AverageExpression(seurat_object_average_expression)
    colnames(avg_expression$SCT)<-paste0("cluster_",colnames(avg_expression$SCT))
    write.table(avg_expression$SCT,gdpath("avrg_expr.tab",root=root,external=root_external),sep="\t",col.names=NA)

    cat("Lets annotate these cells using the gsva database...")
    #lets annotate the cells automatically
    source(gdpath("Masters/project_6_scRNA_seq/analyze_sc_old/cell_type_annotation/obtains_GSVA_for_MatrixColumns_Rvers.R"))
    gsea_scores<-calc_gsva(
      #for InfileMat this weird expression is used because the output of the previous function to get average expression places results in a subfolder
      InfileMat=gdpath("avrg_expr.tab",root=root,external=root_external),#,gdpath(paste0("avrg_expr_gsva","/",project_name,".SEURAT_AverageGeneExpressionPerCluster_cluster_label.tsv"),root=root),
      InfileGmt=infile_gmt,
      PrefixOutfiles=project_name,
      PvalueCutoff=0.10,
      FdrCutoff=0.05,
      output_dir=gdpath("gsva_annotate_cell_types",root=root,external=root_external)
    )

    #get max value of each column which correpsonds to index of the rowname with the cell type
    annotate_cells<-apply(gsea_scores,2,which.max)
    #order clusters numerically (gsva software scarmbles the order)
    annotate_cells<-annotate_cells[order(as.numeric(gsub("cluster_","",names(annotate_cells))))]
    annotate_cells<-rownames(gsea_scores)[annotate_cells]
    annotate_cells<-make.unique(annotate_cells)
    print(annotate_cells)
    return(annotate_cells)

}

read_cluster_mod_data<-function(path_to_mod_data){
	return(read.table(path_to_mod_data,sep="\t",header=TRUE,stringsAsFactors=FALSE))
}

modify_clusters<-function(seurat_object,cluster_mod_data,subroot,root_external){
	if (all(levels(seurat_object$integrated_clusters) %in% cluster_mod_data$original)){
		cat("Original data matches provided in columns... \n")
		#check if there is a manual column
		if ("manual" %in% colnames(cluster_mod_data)){
			cat("Manaul cluster annotation detected...\n")
			cat("Checking manual annotations...\n")
			if (length(cluster_mod_data$manual)== length(cluster_mod_data$original)) {
				#lets remove any clusters with the name "remove"
				seurat_object<-seurat_object[,!(seurat_object$integrated_clusters %in% cluster_mod_data$original[cluster_mod_data$manual=="remove"])] #only keep cells that are not named "remove"
				cluster_mod_data<-cluster_mod_data[!(cluster_mod_data$manual=="remove"),] #in the mod data remove the rows of cells to be removed
				#make sure the remaining cluster_mod_data is correct
				if(all(0:max(as.numeric(cluster_mod_data$manual)) %in% cluster_mod_data$manual)){
					cat("Manual annotation clusters are acceptable. Proceeding with use of manaul annotations.\n")
				}else{
					cat("Error in manual annotation. Ther clusters must be in consecutive integer order. double check. \n")
					quit()
				}
				#modify the seurat clusters to column 2 which contains the manual change
				integrated_clusters<-levels(seurat_object$integrated_clusters)[seurat_object$integrated_clusters] #extract the integrated cluster as charACTER
				for (i in 1:nrow(cluster_mod_data)){
					cat(paste0("converting ",cluster_mod_data$original[i]," to ",cluster_mod_data$manual[i],"\n"))
					integrated_clusters[integrated_clusters==cluster_mod_data$original[i]]<-cluster_mod_data$manual[i]
				}
				names(integrated_clusters)<-rownames(seurat_object@meta.data)
				integrated_clusters<-as.factor(as.numeric(integrated_clusters))
				seurat_object<-AddMetaData(seurat_object,metadata=integrated_clusters,col.name="integrated_clusters")
				Idents(seurat_object)<-seurat_object$integrated_clusters
				#write a log stating that clusters were manually changed
				fileConn<-file(gdpath(paste0("manual_cluster_mod_log.txt"),root=subroot,external=root_external))
			    writeLines(c("Cluster modified by manual file",paste0(cluster_modify_info)), fileConn)
			    close(fileConn)
			    warnings()
				return(seurat_object)
			}else{
				cat("Manual annotation clusters are not correct format. Names must be numbered integers increasing by 1. i.e cluster 1,2,4 are invalid and should be 1,2,3.\n")
				quit()
			}
		}else{
			cat("No manual annotation detected...annotations remain unchanged\n")
			return(seurat_object)
		}
	}else{
		cat("Original data does not match the seurat object. Aborting modification... \n")
		quit()
	}
}

#this function accepts a seurat object containing all cells and a vector of barcodes of recipient cells which to remove. 
remove_cells<-function(seurat_object,barcodes_of_cells_to_be_removed){
	seurat_barcodes<-rownames(seurat_object@meta.data)
	if (!(all(barcodes_of_cells_to_be_removed %in% seurat_barcodes))){
		cat("Warning: one of the recipient barcodes cannot be found in the barcodes of the seurat object. Only barcodes intersecting from both objects will be used to exclude cells\n")
	}
	cells_to_remove<-intersect(seurat_barcodes,barcodes_of_cells_to_be_removed)
	rem<-seurat_object[,!(seurat_barcodes %in% cells_to_remove)]
	return(rem)
}

#read in gavin table
read_in_gavin_table<-function(path){
	return(read.table(file=path,sep=",",header=TRUE,stringsAsFactors=FALSE))
}

# this function removes recipient cells from a seurat object if given gavin_table and the label for the recipient cells used in the 'best_singlet'
# column of the gavin table
remove_recipient_cells<-function(seurat_object,gavin_table,recip_label_in_label_column){
	if (!(any(recip_label_in_label_column %in% gavin_table$label))){
		cat("error. Cannot find at least one instance of recip_label_in_label_column label in the 'best_singlet' column of gavin_table.\n")
		quit()
	}
	recipient_barcodes<-gavin_table$barcode[gavin_table$label==recip_label_in_label_column]
	return(remove_cells(seurat_object=seurat_object,barcodes_of_cells_to_be_removed=recipient_barcodes))
}

#this function adds several meta data for a seurat object pertinant to donor and recipient cell origins
add_deconvolution_info<-function(seurat_object,gavin_table,recip_label_in_label_column){
	cat("conducting donor recipient deconvolution...\n")
	#make sure there is at least one recipient cell
	if (!(any(recip_label_in_label_column %in% gavin_table$label))){
		cat("error. Cannot find at least one instance of recip_label_in_label_column label in the 'best_singlet' column of gavin_table.\n")
		quit()
	}
	#extract barcodes of cells called as the recipient as defined by the recip_label_in_label_column
	recipient_barcodes<-gavin_table$barcode[gavin_table$label==recip_label_in_label_column]
	cat(paste0(length(recipient_barcodes)," recipient barcodes found from Gavin's method. \n"))
	#clean up barcodes
	seurat_object_barcodes<-sapply(rownames(seurat_object@meta.data),FUN=function(x){substr(x,start=1,stop=(str_locate(x,"-")-1))})
	recipient_barcodes_intersecting<-intersect(recipient_barcodes,seurat_object_barcodes)
	#create vector same length as barcodes
	recip_origin<-rep(FALSE,nrow(seurat_object@meta.data))
	#lets replace values in recip_origin with 'TRUE' at locations corresponding to recipient derived cells
	recip_origin[seurat_object_barcodes %in% recipient_barcodes_intersecting & seurat_object@meta.data$timepoint=="post"]<-TRUE
	cat(paste0(sum(recip_origin))," intersecting barcodes identified as recipient.\n")
	#add this vector as a meta data
	seurat_object_modified <- AddMetaData(
	    object = seurat_object,
	    metadata = recip_origin,
	    col.name = 'recipient_origin'
    )
    #add a combined metadata of clusters and recipient origin combined
    seurat_object_modified <- AddMetaData(
	    object = seurat_object_modified,
	    metadata = paste0(seurat_object_modified@meta.data$integrated_clusters,"_",as.character(seurat_object_modified@meta.data$recipient_origin)),
	    col.name = 'integrated_clusters_recipient_origin'
    )
    #add a combined metadata of timepoint and recipient origin combined
    seurat_object_modified <- AddMetaData(
	    object = seurat_object_modified,
	    metadata = paste0(seurat_object_modified@meta.data$timepoint,"_",as.character(seurat_object_modified@meta.data$recipient_origin)),
	    col.name = 'timepoint_recipient_origin'
    )
	return(seurat_object_modified)
}

# this function is similar to "read_in_gavin_table" but we have modified it so it just returns the barcodes of recipient cells
# it was designed for the combined analysis with single cell data from multiple transplants in mind
# this function is meant to be used with "combined_add_deconvolution_info"
extract_recip_barcodes_from_gavin_table<-function(path,label_in_label_column){
	gavin_table<-read.table(file=path,sep=",",header=TRUE,stringsAsFactors=FALSE)
	#make sure there is at least one recipient cell
	if (!(any(label_in_label_column %in% gavin_table$label))){
		cat("error. Cannot find at least one instance of label_in_label_column label in the 'best_singlet' column of gavin_table.\n")
		quit()
	}
	#extract barcodes of cells called as the recipient as defined by the recip_label_in_label_column
	barcodes<-gavin_table$barcode[gavin_table$label==label_in_label_column]
	barcodes<-paste0(barcodes,"-1")
	cat(paste0(length(barcodes),"  barcodes found from Gavin's method. \n"))
	return(barcodes)
}


# it was designed for the combined analysis with single cell data from multiple transplants in mind
#this function adds several meta data for a seurat object pertinant to donor and recipient cell origins
combined_add_deconvolution_info<-function(seurat_object,recipient_barcodes){
	cat("conducting donor recipient deconvolution...\n")
	cat(paste0(length(recipient_barcodes)," total recipient barcodes found from Gavin's method. \n"))
	#clean up barcodes, nevermind maybe we dont have to do anymore
	# seurat_object_barcodes<-sapply(rownames(seurat_object@meta.data),FUN=function(x){substr(x,start=1,stop=(str_locate(x,"-")-1))})
	seurat_object_barcodes<-rownames(seurat_object@meta.data)
	recipient_barcodes_intersecting<-intersect(recipient_barcodes,seurat_object_barcodes)
	#create vector same length as barcodes
	recip_origin<-rep(FALSE,nrow(seurat_object@meta.data))
	#lets replace values in recip_origin with 'TRUE' at locations corresponding to recipient derived cells
	recip_origin[seurat_object_barcodes %in% recipient_barcodes_intersecting & seurat_object@meta.data$timepoint=="post"]<-TRUE
	cat(paste0(sum(recip_origin))," intersecting barcodes (with our droplets call by CellRanger) identified as recipient.\n")
	#add this vector as a meta data
	seurat_object_modified <- AddMetaData(
	    object = seurat_object,
	    metadata = recip_origin,
	    col.name = 'recipient_origin'
    )
    #add a combined metadata of clusters and recipient origin combined
    seurat_object_modified <- AddMetaData(
	    object = seurat_object_modified,
	    metadata = paste0(seurat_object_modified@meta.data$integrated_clusters,"_",as.character(seurat_object_modified@meta.data$recipient_origin)),
	    col.name = 'integrated_clusters_recipient_origin'
    )
    #add a combined metadata of timepoint and recipient origin combined
    seurat_object_modified <- AddMetaData(
	    object = seurat_object_modified,
	    metadata = paste0(seurat_object_modified@meta.data$timepoint,"_",as.character(seurat_object_modified@meta.data$recipient_origin)),
	    col.name = 'timepoint_recipient_origin'
    )
	return(seurat_object_modified)
}

#adds a cell meta data column with the sex of cells. seurat_object$recipient_origin must already be in the object. Also needs "donor_sex" and "recipient sex" in the project data form
add_cell_sex<-function(seurat_object,project_data){
	cell_sex<-rep("blank",nrow(seurat_object@meta.data))
	for (case in unique(seurat_object$ltx_case)){
		print(case)
		#in the cell_sex vector replace the donor cells of the case with the sex (just use first occurence from project data as this is the same whether using pre or post sample)
		cell_sex[seurat_object$ltx_case==case & seurat_object$recipient_origin==FALSE]<-project_data[match(paste0("ltx_",case),project_data$directory),"donor_sex"]
		cell_sex[seurat_object$ltx_case==case & seurat_object$recipient_origin==TRUE]<-project_data[match(paste0("ltx_",case),project_data$directory),"recipient_sex"]
	}
	names(cell_sex)<-rownames(seurat_object@meta.data)
	seurat_object<-AddMetaData(seurat_object,metadata=cell_sex,col.name="cell_sex")
	return(seurat_object)
}