# this function runs the anchored integration from Seurat
# 
# combined_data_list a list of seurat objects to be integrated
# reference_index The index of the element in combined_data_list which should be used as reference for integration
# nfeatures The number of features which should be used for integration.
run_standard_integration<-function(combined_data_list,reference_index,nfeatures=2000){
    combined_data.features <- SelectIntegrationFeatures(object.list = combined_data_list, nfeatures = 2000)
    #un PCA for rpca
    combined_data_list<-lapply(X = combined_data_list, FUN = function(x) {
        # x <- ScaleData(x, features = features, verbose = FALSE) #can remove RunPCA pulls from scale.data slot of SCT
        x <- RunPCA(x, features = combined_data.features, verbose = FALSE)
    })
    combined_data_list <- PrepSCTIntegration(object.list = combined_data_list, anchor.features = combined_data.features, verbose = FALSE)
    combined_data.anchors <- FindIntegrationAnchors(object.list = combined_data_list,reference=reference_index,anchor.features=combined_data.features, normalization.method = "SCT", verbose = TRUE,reduction="rpca")
    combined_data.integrated <- IntegrateData(anchorset = combined_data.anchors, normalization.method = "SCT", 
      verbose = TRUE)
    return(combined_data.integrated)
}