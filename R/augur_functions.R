library("Augur")

#this function runs Augur and writes out a table of results
run_augur_analysis<-function(seurat_object,root,root_external){
	augur = calculate_auc(seurat_object,cell_type_col = "integrated_clusters", label_col = "timepoint",n_threads=10)
	sink(file=gdpath(paste0("cell_type_and_timepoint_distribution_",".txt"),root=root,external=root_external))
	print(table(seurat_object$integrated_clusters,seurat_object$timepoint))
	sink()
	sink(file=gdpath(paste0("augur_results_summary",".txt"),root=root,external=root_external))
	print(augur$AUC)
	sink()
	saveRDS(augur,file=gdpath(paste0("augur_analysis",".rds"),root=root,external=root_external))
}

#augur output some clusters missing because too few samples in the pre timepoint
#$AUC
# # A tibble: 7 x 2
#   cell_type   auc
#   <fct>     <dbl>
# 1 0         0.801
# 2 9         0.745
# 3 8         0.736
# 4 1         0.685
# 5 4         0.675
# 6 6         0.667
# 7 7         0.637

#here we see clusters with very few cells pre are missing from the augur result
# > table(seurat_object$integrated_clusters,seurat_object$timepoint)
    
#      post  pre
#   0  1532 1305
#   1   422  915
#   2   543    4
#   3   516    3
#   4   164  182
#   5   307    8
#   6   200   42
#   7   186  128
#   8   378  147
#   9   106   68
#   10  140   16
#   11  101   14
#   12   98   14