# This script computes differential expression between donor and recipient cells of the same type. It then takes the list of differentially expressed genes and does enrichment analysis. Outputs tables of differentially expressed genes and pathway analysis results to an output folder.

suppressPackageStartupMessages(library(renv))
renv::activate()
suppressPackageStartupMessages(library(optparse))
library(Seurat)
library(trqwe)
source("gdpath.R")
gdlib(c("make_rnk","short_path","GSEA_Pipeline"))
source(gdpath("Masters/project_6_scRNA_seq/scRNA_seq_Cell_Tracker/pathway_analysis/pathway_analysis_functions.R"))
source(gdpath("Masters/project_6_scRNA_seq/scRNA_seq_Cell_Tracker/combined_analysis/combined_resolution_finder.R"))


option_list <- list(
    make_option(c("-p", "--seurat_object_path"), default="NA",
              help="path to seurat object"),
    #
    make_option(c("-a", "--case_list"), default="NA",
              help="a list of cases to use for differential gene expressioon"),
    #
    make_option(c("-b", "--project_data_path"), default="NA",
              help="path to tab delimited file with sample information."),
    #
    make_option(c("-e", "--clusters_with_recip_cells"), default="NA",
              help="cluster numbers containing recip cells seperated by commas"),
    #
    make_option(c("-c", "--geneset_database_pathways"), default="NA",
              help="path to gmt data base for enrichment analysis of pathways"),
    #
    make_option(c("-d", "--geneset_database_tf"), default="NA",
              help="path to gmt data base for enrichment analysis of transcription factors"),
    #
    make_option(c("-f", "--thresh"), default="NA",
              help="percent threshold 0-1 that must be expressed in either ident.1 or ident.2 for gene to be kept in rank list for pathway analysis"),
    #
    make_option(c("-v", "--root"), default="NA",
              help="The directory where the various result directories for each part of the pipeline will be written to."),
    #
    make_option(c("-q", "--root_external"), default="NA",
              help="TRUE or FALSE. If the root is contatined in the cloud synced folder.")
  )

opt <- parse_args(OptionParser(option_list=option_list))


#project_data_path
project_data_path<-opt$project_data_path
cat(paste0("reading in project data at: ",project_data_path," \n"))
project_data<-read.csv(project_data_path,header=TRUE)

#extract arguments
seurat_object_path<-opt$seurat_object_path
root<-opt$root
case_list<-unlist(strsplit(opt$case_list,split=","))
geneset_database_pathways<-opt$geneset_database_pathways
geneset_database_tf<-opt$geneset_database_tf
thresh<-as.numeric(opt$thresh)
recip_containing_cluser<-as.numeric(unlist(strsplit(opt$clusters_with_recip_cells,split=",")))

case_name<-paste0(case_list,collapse="_")

# parameters for writing output
if(opt$root_external=="TRUE"){
  root_external<-TRUE
}else{
  root_external<-FALSE
  cat("internal root detected.\n")
}


#make output directory results directories
dir.create(paste0(root,"deconvolution_analysis/"))
cat("Creating output results directories...\n")
subroot<-paste0(root,"deconvolution_analysis/pathway_analysis/")
case_subroot<-paste0(subroot,paste0(case_list,collapse="_"),"/")
suppressWarnings(dir.create(subroot))
suppressWarnings(dir.create(case_subroot))

cat("reading in seurat object...\n")
seurat_object<-mcreadRDS(seurat_object_path)

#pathway analysis
rank_files<-c()
rank_files_names<-c()
de_result<-list()
de_result_names<-c()
x<-1
for (cell_type in recip_containing_cluser){
	skipped_de<-0
	recip_cells_post<-rownames(seurat_object@meta.data)[(seurat_object@meta.data$integrated_clusters_recipient_origin==paste0(cell_type,"_","TRUE")) & (seurat_object@meta.data$timepoint=="post")]
	donor_cells_post<-rownames(seurat_object@meta.data)[(seurat_object@meta.data$integrated_clusters_recipient_origin==paste0(cell_type,"_","FALSE")) & (seurat_object@meta.data$timepoint=="post")]
	donor_cells_pre<-rownames(seurat_object@meta.data)[(seurat_object@meta.data$integrated_clusters_recipient_origin==paste0(cell_type,"_","FALSE")) & (seurat_object@meta.data$timepoint=="pre")]
	if(file.exists(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_pre.rds"),root=case_subroot,external=root_external)) & file.exists(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_post.rds"),root=case_subroot,external=root_external))){
		cat(paste0("found existing save files. Loading them...\n",cell_type,"\n"))
		compare_recip_to_donor_pre<-mcreadRDS(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_pre.rds"),root=case_subroot,external=root_external))
		compare_recip_to_donor_post<-mcreadRDS(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_post.rds"),root=case_subroot,external=root_external))
	}else{
		cat("running differential gene expression...\n")
		if(length(recip_cells_post)>20){
			#conduct differential expression
			# compare_recip_to_donor_pre<-FindMarkers(seurat_object, ident.1 = recip_cells_post, ident.2 = donor_cells_pre,min.pct=0.10,logfc.threshold=0)
			# compare_recip_to_donor_post<-FindMarkers(seurat_object, ident.1 = recip_cells_post, ident.2 = donor_cells_post,min.pct=0.10,logfc.threshold=0)
			#UNLIKE FINDMARKERS (which does ident.1-ident.2). The more conventional format used, in my function is ident.2-ident.1. THIS IS WHY ORDER
			# IS DIFFERENT BETWEEN ident.1 and ident.2 WHEN USING FINDMARKERS AND MY CUTSOM BPSC function
			# compare_recip_to_donor_pre<-run_BPSC_de(seurat_object=seurat_object,ident.1=donor_cells_pre,ident.2=recip_cells_post,pseudocount.use=1)
			# compare_recip_to_donor_post<-run_BPSC_de(seurat_object=seurat_object,ident.1=donor_cells_post,ident.2=recip_cells_post,pseudocount.use=1)
			# compare_recip_to_donor_pre<-run_de_limma_voom(seurat_object=seurat_object,ident.1=donor_cells_pre,ident.2=recip_cells_post)
			# compare_recip_to_donor_post<-run_de_limma_voom(seurat_object=seurat_object,ident.1=donor_cells_post,ident.2=recip_cells_post)
			# compare_recip_to_donor_pre<-FindMarkers(seurat_object,ident.1=donor_cells_pre,ident.2=recip_cells_post,test.use="negbinom",latent.vars = "ltx_case",min.pct = 0.10,logfc.threshold=0)
			# compare_recip_to_donor_post<-FindMarkers(seurat_object,ident.1=donor_cells_post,ident.2=recip_cells_post,test.use="negbinom",latent.vars = "ltx_case",min.pct = 0.10,logfc.threshold=0)
			compare_recip_to_donor_pre<-run_nebula_raw_counts_deconvolution(seurat_object=seurat_object,ident.1=donor_cells_pre, ident.2=recip_cells_post,metadata=metadata)
			compare_recip_to_donor_post<-run_nebula_raw_counts_deconvolution(seurat_object=seurat_object,ident.1=donor_cells_post, ident.2=recip_cells_post,metadata=metadata)
			mcsaveRDS(compare_recip_to_donor_pre,gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_pre.rds"),root=case_subroot,external=root_external))
			mcsaveRDS(compare_recip_to_donor_post,gdpath(paste0("cluster_",cell_type,"_",case_name,"_","compare_recip_to_donor_post.rds"),root=case_subroot,external=root_external))
		}else{
			cat(paste0("\n\n###################################################################\n"))
			cat(paste0(cell_type," did not have at least 50 cells in the recipient group.\n"))
			cat(paste0("###################################################################\n"))
			skipped_de<-1
		}
	}
	if (skipped_de!=1){
		de_result[[x]]<-compare_recip_to_donor_pre
		de_result[[x+1]]<-compare_recip_to_donor_post
		x<-x+2
		de_result_names<-c(de_result_names,paste0(cell_type,c("_recip_to_donor_pre","_recip_to_donor_post")))
		# #run_BPSC_de outputs a different format
		# compare_recip_to_donor_pre_table<-data.frame(compare_recip_to_donor_pre$PVAL,compare_recip_to_donor_pre$logFC)
		# colnames(compare_recip_to_donor_pre_table)<-c("PVAL","logFC")
		# compare_recip_to_donor_post_table<-data.frame(compare_recip_to_donor_post$PVAL,compare_recip_to_donor_post$logFC)
		# colnames(compare_recip_to_donor_post_table)<-c("PVAL","logFC")
		# #remove genes that dont meet our threshold for minimum expression of at least 1 copy in x percent of cells in either ident.1 or ident.2
		# compare_recip_to_donor_pre_table<-compare_recip_to_donor_pre_table[compare_recip_to_donor_pre$ident.1.pct>thresh | compare_recip_to_donor_pre$ident.2.pct>thresh ,]
		# compare_recip_to_donor_post_table<-compare_recip_to_donor_post_table[compare_recip_to_donor_post$ident.1.pct>thresh | compare_recip_to_donor_post$ident.2.pct>thresh,]
		#lets get fdr corrected values
		# compare_recip_to_donor_pre_table$fdr<-p.adjust(compare_recip_to_donor_pre_table[,"PVAL"],method="fdr")
		# compare_recip_to_donor_post_table$fdr<-p.adjust(compare_recip_to_donor_post_table[,"PVAL"],method="fdr")
		#make histogram plots
		cat("making histogram plots\n")
		pdf(gdpath(paste0("cluster_",cell_type,"_",case_name,"_","pvalue_distribution.pdf"),root=case_subroot,external=root_external))
		print(hist(compare_recip_to_donor_pre[,"p_val"],breaks=20,main=paste0(cell_type,": recipient_post vs donor_pre")))
		print(hist(compare_recip_to_donor_post[,"p_val"],breaks=20,main=paste0(cell_type,": recipient_post vs donor_post")))
		dev.off()

		#write tab files
		write.table(compare_recip_to_donor_pre,file=gdpath(paste0(cell_type,"_recipient_post_vs_donor_pre_DE",".tab"),root=case_subroot,external=root_external),sep="\t",row.names=TRUE,col.names=NA)
		write.table(compare_recip_to_donor_post,file=gdpath(paste0(cell_type,"_recipient_post_vs_donor_post_DE",".tab"),root=case_subroot,external=root_external),sep="\t",row.names=TRUE,col.names=NA)	

		#make sure only keep genes that meet thresh
		compare_recip_to_donor_pre<-compare_recip_to_donor_pre[(compare_recip_to_donor_pre[,"ident.1.pct"]>thresh | compare_recip_to_donor_pre[,"ident.2.pct"]>thresh),]
		compare_recip_to_donor_post<-compare_recip_to_donor_post[(compare_recip_to_donor_post[,"ident.1.pct"]>thresh | compare_recip_to_donor_post[,"ident.2.pct"]>thresh),]

		#make rnk files
		cat("making rnk files\n")
		# make_rnk_simple_by_fold_change(compare_recip_to_donor_pre[,c("p_val","avg_logFC")],gdpath(paste0(cell_type,"_recipient_post_vs_donor_pre",".rnk"),root=case_subroot,external=root_external),clean=F)
		# make_rnk_simple_by_fold_change(compare_recip_to_donor_post[,c("p_val","avg_logFC")],gdpath(paste0(cell_type,"_recipient_post_vs_donor_post",".rnk"),root=case_subroot,external=root_external),clean=F)
		make_rnk_simple(compare_recip_to_donor_pre[,c("p_val","avg_logFC")],gdpath(paste0(cell_type,"_recipient_post_vs_donor_pre",".rnk"),root=case_subroot,external=root_external),clean=F)
		make_rnk_simple(compare_recip_to_donor_post[,c("p_val","avg_logFC")],gdpath(paste0(cell_type,"_recipient_post_vs_donor_post",".rnk"),root=case_subroot,external=root_external),clean=F)
		rank_files<-append(rank_files,gdpath(paste0(cell_type,"_recipient_post_vs_donor_pre",".rnk"),root=case_subroot,external=root_external))
		rank_files<-append(rank_files,gdpath(paste0(cell_type,"_recipient_post_vs_donor_post",".rnk"),root=case_subroot,external=root_external))
		rank_files_names<-append(rank_files_names,paste0(cell_type,"_recipient_post_vs_donor_pre"))
		rank_files_names<-append(rank_files_names,paste0(cell_type,"_recipient_post_vs_donor_post"))
		# write.table(compare_recip_to_donor_pre,file=gdpath(paste0(cell_type,"_recipient_post_vs_donor_pre_DE",".tab"),root=subroot,external=root_external),sep="\t",row.names=TRUE,col.names=NA)
		# write.table(compare_recip_to_donor_post,file=gdpath(paste0(cell_type,"_recipient_post_vs_donor_post_DE",".tab"),root=subroot,external=root_external),sep="\t",row.names=TRUE,col.names=NA)
	}
}

names(de_result)<-de_result_names
mcsaveRDS(de_result,gdpath(paste0("all_comparisons_recip_to_donor_post.rds"),root=case_subroot,external=root_external))

names(rank_files)<-rank_files_names


#run pathway analysis comparing donor and recipient cells with a pathway database
run_GSEA(
    rank_files=rank_files,
    gmt_database=geneset_database_pathways,
    out_folder=case_subroot,
    root=root,
    root_external=root_external,
    project_name=project_name,
    min_geneset=8,
    max_geneset=300,
    output_append="pathways"
)

#run pathway analysis comparing donor and recipient cells with a transcription factor database
run_GSEA(
    rank_files=rank_files,
    gmt_database=geneset_database_tf,
    out_folder=case_subroot,
    root=root,
    root_external=root_external,
    project_name=project_name,
    min_geneset=8,
    max_geneset=300,
    output_append="TF"
)