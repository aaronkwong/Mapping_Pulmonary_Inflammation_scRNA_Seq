# takes a seurat object (which should contain only a single cluster of cells) and increases the clustering resolution untill the desired number of subclusters are found. Returns the clustering resolution that was used. 
# @param seurat_object A seurat object containing cells to be subclustered
# @param desired_clusters The number of sub populations to find
# @param n_pcs Number of principal components to use when sub clustering
find_clusters_by_communities<-function(seurat_object,desired_clusters,n_pcs){
	seurat_object<-RunPCA(seurat_object,verbose=FALSE)
	seurat_object<-FindNeighbors(seurat_object,dims=1:n_pcs,verbose=FALSE)
	found<-0
	tested_resolutions<-c(0,0.1)
	while(TRUE){
		Sys.sleep(0.1)
		seurat_object<-FindClusters(seurat_object,resolution=tested_resolutions[2],verbose=FALSE)
		clusters_made<-length(unique(seurat_object@meta.data[,paste0("integrated_snn_res.",tested_resolutions[2])]))
		print(clusters_made)
		if(clusters_made==desired_clusters){
			cat(paste0("found desired clusters\n"))
			return(tested_resolutions[2])
		}else if(clusters_made<desired_clusters){
			if(tested_resolutions[1]>tested_resolutions[2]){
				tested_resolutions<-c(tested_resolutions,(tested_resolutions[1]+tested_resolutions[2])/2)[-1]
			}else{
				tested_resolutions<-c(tested_resolutions,tested_resolutions[2]*1.2)[-1]
			}
		}else if(clusters_made>desired_clusters){
			if(tested_resolutions[1]>tested_resolutions[2]){
				tested_resolutions[1]<-0
			}else{
				tested_resolutions<-c(tested_resolutions,(tested_resolutions[1]+tested_resolutions[2])/2)[-1]
			}
		}
	}
}

# 

# @param anchored_data seurat object containing anchored integrated cells
# @param cluster_for_subclustering the name (Idents) of the cluster to be subclustered
# @param desired_subclusters The number of sub communities you want to find from the cluster
# @param n_pcs number of PCs to use
subcluster<-function(anchored_data,cluster_for_subclustering,desired_subclusters,n_pcs){
	cat(paste0("splitting cluster ",cluster_for_subclustering," into ",desired_subclusters," clusters.\n"))
	extracted_cluster<-anchored_data[,anchored_data$integrated_clusters==cluster_for_subclustering]
	#using the appropriate resolution cluster
	required_resolution<-find_clusters_by_communities(seurat_object=extracted_cluster,desired_clusters=desired_subclusters,n_pcs=n_pcs)
	extracted_cluster<-RunPCA(extracted_cluster,verbose=FALSE)
	extracted_cluster<-FindNeighbors(extracted_cluster,dims=1:n_pcs,verbose=FALSE)
	extracted_cluster<-FindClusters(extracted_cluster,resolution=required_resolution,verbose=FALSE)

	for (subcluster in 1:(desired_subclusters-1)){
		new_cluster_number<-max(as.numeric(levels(anchored_data$integrated_clusters)))+1
		#get barcodes of subcluster
		barcodes_new_cluster<-rownames(extracted_cluster@meta.data)[extracted_cluster@meta.data[,paste0("integrated_snn_res.",required_resolution)]==subcluster] # the other half labelled cluster 0 will keep the same cluster number
		# barcodes_new_cluster<-rownames(extracted_cluster@meta.data)[extracted_cluster@meta.data$seurat_clusters==subcluster] # the other half labelled cluster 0 will keep the same cluster number
		cat(paste0(length(barcodes_new_cluster)," cells subclustered out.\n"))
		#extract all barcodes from anchored object
		integrated_clusters<-as.numeric(levels(anchored_data$integrated_clusters)[anchored_data$integrated_clusters])
		cat(paste0("matching cells ",sum(rownames(anchored_data@meta.data) %in% barcodes_new_cluster),"\n"))
		#using barcodes of the new subcluster, find these cells in the anchored object and reassign the clustering
		integrated_clusters[rownames(anchored_data@meta.data) %in% barcodes_new_cluster]<-new_cluster_number
		names(integrated_clusters)<-rownames(anchored_data@meta.data)
		integrated_clusters<-as.factor(integrated_clusters)
		levels(integrated_clusters)<-as.character(seq(from=0,to=length(unique(integrated_clusters))-1,by=1))
		#reappend the modified cluster data to the original anchored object. Repeat as necessary if desired subclusters >2
		anchored_data<-AddMetaData(anchored_data,metadata = integrated_clusters,col.name="integrated_clusters")
	}
	return(anchored_data)
}
