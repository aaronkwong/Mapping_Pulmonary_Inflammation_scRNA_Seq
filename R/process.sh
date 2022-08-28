# conduct the first simple integration
Rscript "${script_dir}driver_scripts/combine_simple.R" \
--root=$root_directory \
--root_external=$root_external \
--pc_use=$n_pcs \
--project_name=$project_name \
--project_data_path="$project_data_path" \
--cluster_resolution=$cluster_resolution_pre_cor \
--project_directory=$project_directory \
--max_cell_scaling=$max_cell_scaling

# conduct the first anchored integration to make clusters for QC
Rscript "${script_dir}driver_scripts/combine_anchored.R" \
--root=$root_directory \
--root_external=$root_external \
--pc_use=$n_pcs \
--project_name=$project_name \
--project_data_path="$project_data_path" \
--cluster_resolution=$cluster_resolution_pre_cor \
--project_directory=$project_directory \
--max_cell_scaling=$max_cell_scaling

# conduct manual cell type annotation-------------------------------------------------------------------------------------------
# takes as input the path to the seurat object (cleaned by the doublet removal above) saved as rds object
# outputs results into a pathway analysis folder with a file with cluster differential expression [gdpath(paste0(project_name,"_all_cluster_DE.rdata"),root=subroot,external=root_external)]
Rscript "${script_dir}driver_scripts/annotate_clusters.R" \
--min_prct_cells=0.2 \
--logfc_threshold=0.25 \
--infile_gmt="$cell_annotation_gmt_path" \
--project_name=$project_name \
--run_deconvolution=$run_deconvolution \
--root=$root_directory \
--root_external=$root_external \
--dir_append=pre_cor \
--project_data_path="$project_data_path" \
--seurat_object_path="initial_integration/combined_data.rds" \
--make_save_rds=FALSE \
--use_manual_anno=FALSE \
--make_marker_table=FALSE

#run background removal from the object generated after the inital integration
Rscript "${script_dir}driver_scripts/deconvolve_background_RNA.R" \
--root=$root_directory \
--root_external=$root_external \
--pc_use=$n_pcs \
--project_name=$combined_analysis

#run doublet removal from the object generated after the inital integration
Rscript "${script_dir}driver_scripts/multiplet_removal.R" \
--root=$root_directory \
--root_external=$root_external \
--pc_use=$n_pcs \
--project_name=$combined_analysis \
--project_data_path="$project_data_path"

#post cor integration pt.1
Rscript "${script_dir}driver_scripts/post_qc_combine_simple.R" \
--root=$root_directory \
--root_external=$root_external \
--pc_use=$n_pcs \
--project_name=$project_name \
--project_data_path="$project_data_path" \
--cluster_resolution=$cluster_resolution_post_cor \
--remove_doublets="TRUE"

#post cor integration pt.2
Rscript "${script_dir}driver_scripts/post_qc_combine_anchored.R" \
--root=$root_directory \
--root_external=$root_external \
--pc_use=$n_pcs \
--project_name=$project_name \
--project_data_path="$project_data_path" \
--cluster_resolution=$cluster_resolution_post_cor \
--remove_doublets="TRUE"

#lets make a quick annotation to help us identify any clusters needing subclustering
Rscript "${script_dir}driver_scripts/annotate_clusters.R" \
--min_prct_cells=0.1 \
--logfc_threshold=0.25 \
--infile_gmt="$cell_annotation_gmt_path" \
--project_name=$project_name \
--run_deconvolution=$run_deconvolution \
--root=$root_directory \
--root_external=$root_external \
--dir_append=subcluster_inspection \
--project_data_path="$project_data_path" \
--seurat_object_path="post_QC_processing/"$project_name"_ANCHORED_integration.rds" \
--make_save_rds=TRUE \
--use_manual_anno=FALSE \
--make_marker_table=TRUE

# now lets do subclustering on any clusters requiring it
Rscript "${script_dir}driver_scripts/subcluster.R" \
--pc_use=$n_pcs \
--project_name=$project_name \
--input_dir="post_QC_processing/" \
--root=$root_directory \
--root_external=$root_external \
--make_save=TRUE

# INTEGRATED AUTO
Rscript "${script_dir}driver_scripts/annotate_clusters.R" \
--min_prct_cells=0.1 \
--logfc_threshold=0.25 \
--infile_gmt="$cell_annotation_gmt_path" \
--project_name=$project_name \
--run_deconvolution=$run_deconvolution \
--root=$root_directory \
--root_external=$root_external \
--dir_append=post_cor_integrated_auto \
--project_data_path="$project_data_path" \
--seurat_object_path="subclustering/project_data_ANCHORED_integration.rds" \
--make_save_rds=TRUE \
--use_manual_anno=FALSE \
--make_marker_table=FALSE


#run annotate cells but dont use any manual annotations. Use this version to merge clusters
# COMBINED AUTO
Rscript "${script_dir}driver_scripts/annotate_clusters.R" \
--min_prct_cells=0.1 \
--logfc_threshold=0.25 \
--infile_gmt="$cell_annotation_gmt_path" \
--project_name=$project_name \
--run_deconvolution=$run_deconvolution \
--root=$root_directory \
--root_external=$root_external \
--dir_append=post_cor_combined_auto \
--project_data_path="$project_data_path" \
--seurat_object_path="subclustering/combined_data.rds" \
--make_save_rds=TRUE \
--use_manual_anno=FALSE \
--make_marker_table=TRUE

#run singleR. It would be ideal to run this after the auto run of integrated and combined. SingleR data can help us justify merging clusters
Rscript "${script_dir}driver_scripts/singleR_annotation.R" \
--root=$root_directory \
--root_external=$root_external \
--travag_reference_path=$singleR_reference_dataset \
--test_cells_path=$root_directory"annotate_cells_post_cor_combined_auto/"$project_name"_man_anno_WITH_recipient_cells.rds"


# INTEGRATED MANUAL
Rscript "${script_dir}driver_scripts/annotate_clusters.R" \
--min_prct_cells=0.1 \
--logfc_threshold=0.25 \
--infile_gmt="$cell_annotation_gmt_path" \
--project_name=$project_name \
--run_deconvolution=$run_deconvolution \
--root=$root_directory \
--root_external=$root_external \
--dir_append=post_cor_integrated \
--project_data_path="$project_data_path" \
--seurat_object_path="subclustering/project_data_ANCHORED_integration.rds" \
--make_save_rds=TRUE \
--use_manual_anno=TRUE \
--make_marker_table=FALSE

# COMBINED MANUAL
Rscript "${script_dir}driver_scripts/annotate_clusters.R" \
--min_prct_cells=0.1 \
--logfc_threshold=0.25 \
--infile_gmt="$cell_annotation_gmt_path" \
--project_name=$project_name \
--run_deconvolution=$run_deconvolution \
--root=$root_directory \
--root_external=$root_external \
--dir_append=post_cor_combined \
--project_data_path="$project_data_path" \
--seurat_object_path="subclustering/combined_data.rds" \
--make_save_rds=TRUE \
--use_manual_anno=TRUE \
--make_marker_table=TRUE

#pathway analysis
Rscript "${script_dir}driver_scripts/DGE_and_pathway_analysis.R" \
--seurat_object_path=$root_directory"annotate_cells_post_cor_combined/"$project_name"_man_anno.rds" \
--min_geneset=10 \
--max_geneset=500 \
--root=$root_directory \
--root_external=$root_external \
--project_name=$project_name \
--geneset_database_pathways_1="${onedrive}Masters/Reference Data/Pathway Analysis/Bader_genesets/Human/Pathways/Human_GOBP_AllPathways_no_GO_iea_March_01_2021_symbol.gmt" \
--geneset_database_pathways_2="NA" \
--geneset_database_pathways_3="NA" \
--geneset_database_tf_1="NA" \
--geneset_database_tf_2="NA" \
--geneset_database_tf_3="${onedrive}Masters/Reference Data/Pathway Analysis/ChEA3/ENCODE_ChIP_seq.gmt" \
--thresh=0.2

#######################################################################################################
# Deconvolution analysis --------------- START --------------------------------------
#######################################################################################################
#visualze deconvoluted cells
Rscript "${script_dir}driver_scripts/deconvolve_cells.R" \
--seurat_object_path=$root_directory"annotate_cells_post_cor_combined/"$project_name"_man_anno_WITH_recipient_cells.rds" \
--clusters_with_recip_cells=$recip_containing_clusters \
--root=$root_directory \
--root_external=$root_external

#run deconvolution analysis all cases combined, cell types 6 and 11
Rscript "${script_dir}driver_scripts/DGE_pathway_analysis_donor_versus_recipient.R" \
--seurat_object_path=$root_directory"annotate_cells_post_cor_combined/"$project_name"_man_anno_WITH_recipient_cells.rds" \
--project_data_path="$project_data_path" \
--case_list="case_2,case_3,case_4,case_5,case_6,case_7" \
--clusters_with_recip_cells=$recip_containing_clusters \
--geneset_database_pathways="${onedrive}Masters/Reference Data/Pathway Analysis/Bader_genesets/Human/Pathways/Human_GOBP_AllPathways_no_GO_iea_March_01_2021_symbol.gmt" \
--geneset_database_tf="${onedrive}Masters/Reference Data/Pathway Analysis/ChEA3/ENCODE_ChIP_seq.gmt" \
--thresh=0.2 \
--root=$root_directory \
--root_external=$root_external

#######################################################################################################
# Deconvolution analysis --------------- END --------------------------------------
#######################################################################################################

#cellchat
Rscript "${script_dir}driver_scripts/cell_interaction_analysis.R" \
--seurat_object_path_with_recip_cells=$root_directory"annotate_cells_post_cor_combined/"$project_name"_man_anno_WITH_recipient_cells.rds" \
--project_name=$project_name \
--root=$root_directory \
--root_external=$root_external \
--clusters_containing_recip_cells=$recip_containing_clusters


#lets run augur
Rscript "${script_dir}driver_scripts/cell_prioritization_driver.R" \
--seurat_object_path=$root_directory"annotate_cells_post_cor_combined/"$project_name"_man_anno.rds" \
--root=$root_directory \
--root_external=$root_external 