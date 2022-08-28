#This script computes differential gene expression for each donor cell type before and after reperfusion using NEBULA ("NEBULA is a fast negative binomial mixed model for differential or co-expression analysis of large-scale multi-subject single-cell data"). The list of differentially expressed genes is then used for enrichment analysis (GSEA). Up to 6 different GMT files can be specified to run multiple enrichment analysis, results from each analysis will be output in their own subdirectory. 

suppressPackageStartupMessages(library(renv))
renv::activate()
suppressPackageStartupMessages(library(optparse))
library(Seurat)
source(gdpath("Masters/project_6_scRNA_seq/scRNA_seq_Cell_Tracker/combined_analysis/combined_resolution_finder.R"))


#source custom functions
source("gdpath.R")
gdlib(c("path_manipulation"))


option_list <- list(
    make_option(c("-r", "--seurat_object_path"), default="NA",
        help="Clustered seurat object, each cluster will have differential expression tested between timepoints"),
    #
    make_option(c("-x", "--min_geneset"), default="NA",
        help="The smallest gene sets to be used for GSEA."),
    #
    make_option(c("-y", "--max_geneset"), default="NA",
        help="The largest gene sets to be used for GSEA."),
    #
    make_option(c("-t", "--root"), default="NA",
        help="The directory where the various result directories for each part of the pipeline will be written to."),  
    #
    make_option(c("-p", "--root_external"), default="NA",
        help="TRUE or FALSE. If the root is contatined in the cloud synced folder."),
    #
    make_option(c("-c", "--project_name"), default="NA",
        help="A name for the project."),
    #
    make_option(c("-a", "--geneset_database_pathways_1"), default="NA",
        help="path to a geneset database (.gmt)"),
    #
    make_option(c("-b", "--geneset_database_pathways_2"), default="NA",
        help="path to a geneset database (.gmt)"),
    #
    make_option(c("-d", "--geneset_database_pathways_3"), default="NA",
        help="path to a geneset database (.gmt)"),
    #
    make_option(c("-e", "--geneset_database_tf_1"), default="NA",
        help="path to a geneset database (.gmt)"),
    #
    make_option(c("-f", "--geneset_database_tf_2"), default="NA",
        help="path to a geneset database (.gmt)"),
    #
    make_option(c("-g", "--geneset_database_tf_3"), default="NA",
        help="path to a geneset database (.gmt)"),
    #
    make_option(c("-i", "--thresh"), default="NA",
        help="Minimum percent expression cutoff for genes to be used in pathway analysis"),
    #
    make_option(c("-z", "--logFC_thresh"), default="NA",
        help="Minimum logFC cutoff for genes to be used in pathway analysis")
  )

# set run params
opt <- parse_args(OptionParser(option_list=option_list))
options(future.globals.maxSize = 120000 * 1024^2)
library(doParallel)
registerDoParallel(cores=8)

#set params
seurat_object_path<-opt$seurat_object_path
max_geneset<-as.numeric(opt$max_geneset)
min_geneset<-as.numeric(opt$min_geneset)
root<-opt$root
geneset_database_pathways<-c(opt$geneset_database_pathways_1,opt$geneset_database_pathways_2,opt$geneset_database_pathways_3)
geneset_database_tf<-c(opt$geneset_database_tf_1,opt$geneset_database_tf_2,opt$geneset_database_tf_3)
thresh<-as.numeric(opt$thresh)
logFC_thresh<-as.numeric(opt$logFC_thresh)

#clean list to only contain databases
geneset_database_pathways<-geneset_database_pathways[geneset_database_pathways!="NA"]
geneset_database_tf<-geneset_database_tf[geneset_database_tf!="NA"]

#print out databases that will be used
cat("Pathway analysis will be run on:\n")
print(geneset_database_pathways)
cat("TF analysis will be run on:\n")
print(geneset_database_tf)

#set write output parameters
if (opt$root_external=="TRUE"){
    root_external<-TRUE
}else{
    root_external<-FALSE
}
project_name<-opt$project_name


if(!(
  # here we only need greater than 5 arguments, because we need at least one database.
    (sum(opt!="NA")>5+1) &
    file.exists(seurat_object_path)
    ))
    {
    cat("driver check failed... One of the arguments is missing, or the RDS file specified does not exists.\n")
    cat(paste0("number of args: ",length(opt),"\n"))
    cat(paste0("seurat_object_exists: ",file.exists(seurat_object_path),"\n"))
    quit()
}

cat("driver check passed....\n")
source(gdpath("Masters/project_6_scRNA_seq/scRNA_seq_Cell_Tracker/pathway_analysis/pathway_analysis_functions.R"))

subroot<-paste0(root,"enrichment_analysis/")

#if differential gene expression was already run, load the list of genes instead of re-computing
if(file.exists(gdpath("rnk_file_names.rds",root=subroot,external=root_external))){
    cat("\n#########################################################\n")
    cat("Past differential gene expression file found. Using it.\n\n")
    #load differentially expressed genes computed for each cluster
    rank_files<-readRDS(gdpath("rnk_file_names.rds",root=subroot,external=root_external))
}else{
    rank_files<-run_pathway_analysis_combined( 
        seurat_object=readRDS(seurat_object_path),
        subroot=subroot,
        root=root,
        root_external=root_external,
        project_name=project_name,
        thresh=thresh
    )

    saveRDS(rank_files,file=gdpath("rnk_file_names.rds",root=subroot,external=root_external)) # saves all genes, no cutoff so that cutoff can change in subsequent runs
}

#lets cutoff each rnk file by a min percentage expressed in either ident.1 or ident.2
rank_files<-lapply(rank_files,FUN=function(x){
        if (x!="insufficient_cells"){
            return(x[(x$ident.1.pct>thresh | x$ident.2.pct>thresh),])
        }else{
            return(x)
        }
    }
)

#run pathway enrichment analysis for each input database
if(length(geneset_database_pathways)!=0){
    for (database in geneset_database_pathways){
      #lets create a subfolder for pathway genesets 
      pathways_subroot<-paste0(subroot,paste0(last_obj_from_path(database),"/"))
      run_GSEA(
          rank_files=rank_files,
          min_geneset=min_geneset,
          max_geneset=max_geneset,
          gmt_database=database,
          out_folder=pathways_subroot,
          root=root,
          root_external=root_external,
          project_name=project_name
      )
    }
}

#run Transcription Factor enrichment analysis
#lets create a subfolder for TF genesets 
if (length(geneset_database_tf)!=0){
    for (database in geneset_database_tf){
      TF_subroot<-paste0(subroot,paste0(last_obj_from_path(database),"/"))
      run_GSEA(
        rank_files=rank_files,
        min_geneset=min_geneset,
        max_geneset=max_geneset,
        gmt_database=database,
        out_folder=TF_subroot,
        root=root,
        root_external=root_external,
        project_name=project_name
        )
    }
}

