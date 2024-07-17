#set run params

results_folder="/mnt/c/phd/project_6_human_ltx/"
#variables that need to change for each patient
root_directory=$results_folder"testing/combined_analysis_v6_rerun_2024_06_01/"
#set a project name
project_name=combined_analysis_v6_working
root_external=TRUE

###################################################################

Rscript "./setup_renv.R"

#pathway analysis
Rscript "./R_Scripts/combined_analysis/driver_pathway_analysis_combined.R" \
--seurat_object_path="./annotate_cells_post_cor_combined/all_cells.rds" \
--cluster_resolution="low" \
--include_recip="FALSE" \
--min_geneset=10 \
--max_geneset=500 \
--root=$root_directory \
--root_external=$root_external \
--project_name=$project_name

#donor recipient differential gene expression using low resolution clusters
Rscript "./R_Scripts/combined_analysis/deconvolution_analysis.R" \
--seurat_object_path="./annotate_cells_post_cor_combined/all_cells.rds" \
--cluster_resolution="low" \
--include_recip="TRUE" \
--case_list="case_2,case_3,case_4,case_5,case_6,case_7" \
--clusters_with_recip_cells="5,6,10,11,17" \
--root=$root_directory \
--root_external=$root_external \
--analysis_name="lowRes"

#donor recipient differential gene expression using high resolution clusters
#run deconvolution analysis all cases combined, cell types 6 and 11
Rscript "./R_Scripts/combined_analysis/deconvolution_analysis.R" \
--seurat_object_path="./annotate_cells_post_cor_combined/all_cells.rds" \
--cluster_resolution="high" \
--include_recip="TRUE" \
--case_list="case_2,case_3,case_4,case_5,case_6,case_7" \
--clusters_with_recip_cells="5,6,10,11,17,18,19" \
--root=$root_directory \
--root_external=$root_external \
--analysis_name="highRes"


#cellchat
Rscript "./R_Scripts/combined_analysis/combined_cellchat.R" \
--seurat_object_path_with_recip_cells"./annotate_cells_post_cor_combined/all_cells.rds" \
--cluster_resolution="low" \
--include_recip="TRUE" \
--project_name=$project_name \
--root=$root_directory \
--root_external=$root_external \
--clusters_containing_recip_cells=$recip_containing_clusters
