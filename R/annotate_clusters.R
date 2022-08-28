#this script takes a seurat object and clusters cells according to specified parameters. Can apply manually specified cluster annotations, or can aumotaically annotate cells if provided with a GMT file containing cell type markers. This script can also incorporate donor and recipient deconvolution data from scSNV. This script will also generate markers for each cluster and write this out as a tab-delimited file. The clustered Seurat object and summary plots are generated.

suppressPackageStartupMessages(library(renv))
renv::activate()
suppressPackageStartupMessages(library(trqwe))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-a", "--min_prct_cells"), default="NA",
        help="The minimum prct of cells cutoff when testing for differential gene expression of cell markers"),
    #
    make_option(c("-b", "--logfc_threshold"), default="NA",
        help="The logFC threshold for testing differential gene expression of cell markers"),
    #
    make_option(c("-c", "--project_name"), default="NA",
        help="A name for the project."),
    #
    make_option(c("-d", "--run_deconvolution"), default="NA",
        help="TRUE or FALSE should scSNV deconvolution data be added?"),
    #
    make_option(c("-f", "--root"), default="NA",
        help="The directory where the various result directories for each part of the pipeline will be written to."),  
    #
    make_option(c("-g", "--root_external"), default="NA",
        help="TRUE or FALSE if the root is contatined in the cloud synced folder."),  
    #
    make_option(c("-i", "--infile_gmt"), default="NA",
        help="GMT file containing annotations for each cell type, used for auto annotation"),  
    #
    make_option(c("-j", "--dir_append"), default="NA",
        help="string to append to the output directory which will be written to the root directory"),
    #
    make_option(c("-k", "--project_data_path"), default="NA",
        help="path to excel file containing data for each sample"),
    #
    make_option(c("-y", "--seurat_object_path"), default="NA",
        help="path to quality controlled data"),  
    #
    make_option(c("-m", "--make_save_rds"), default=FALSE,
        help="TRUE or FALSE. Should saves of the output be made?"),  
    #
    make_option(c("-l", "--use_manual_anno"), default=FALSE,
        help="TRUE or FALSE. If TRUE, manual annotations (as specified in manual_changes/manual_cluster_names.tab) will be used for annotation. If FALSE, auto annotation will be used."),  
    #
    make_option(c("-n", "--make_marker_table"), default=FALSE,
        help="TRUE or FALSE. Should markers for each cell cluster be computed and output?")
    )

opt <- parse_args(OptionParser(option_list=option_list))

#source custom functions
source("gdpath.R")

#set resource allocation
options(future.globals.maxSize = 60000 * 1024^2)
plan("multiprocess", workers = 16)
set.seed(42)

#set analysis parameters
#minimum percent cutoff when conducting differential gene expression to find cluster markers 
min_prct_cells<-as.numeric(opt$min_prct_cells)
#minimum logFC threshold when conducting differential gene expression to find cluster markers 
logfc.threshold<-as.numeric(opt$logfc_threshold)
#logical whether to run donor recipient deconvolution data. convert string to logical
if (opt$run_deconvolution=="TRUE"){
    run_deconvolution<-TRUE
}else if(opt$run_deconvolution=="FALSE"){
    run_deconvolution<-FALSE
}else{
    cat("run deconvolution is neither 'TRUE' nor 'FALSE'. \n")
    quit()
}
#path to .tab file containg input information for deconvolution
deconvolution_table_path<-opt$deconvolution_table_path
#root directory where results are written
root<-opt$root
#project name
project_name<-opt$project_name
#logical whether root is in a cloud synced folder or not. convert string to logical
if (opt$root_external=="TRUE"){
    root_external<-TRUE
}else{
    root_external<-FALSE
}
#gmt file containing cell type and the gene set markers used to annotate each cell ty
infile_gmt<-opt$infile_gmt
#prefix to be appended to results as this script is run at various stages of the pipleline
dir_append<-opt$dir_append
#path to .tab file containg project data
project_data_path<-opt$project_data_path
#the seurate object to analyze
seurat_object_path<-opt$seurat_object_path
if (opt$make_save_rds=="TRUE"){
    cat("Save objects will be created.\n")
  make_save_rds<-TRUE
}else{
    make_save_rds<-FALSE
}

#logical whether to use manual annotations as defined in paste0(root,"/manual_changes/manual_cluster_names.tab")
if (opt$use_manual_anno=="TRUE"){
    cat("manual anno will be used.\n")
    use_manual_anno<-TRUE
}else{
    cat("manual anno not used.\n")
    use_manual_anno<-FALSE
}

#logical where to compute cluster markers and output a table
if (opt$make_marker_table=="TRUE"){
    cat("marker tables will be created.\n")
    make_marker_table<-TRUE
}else{
    cat("marker tables will not be created.\n")
    make_marker_table<-FALSE
}

#logical if annotated seurat_object should be saved as a rds object 
if(!make_save_rds==FALSE){
    if(make_save_rds=="TRUE"){
        make_save_rds<-TRUE
    }else{
        cat("make_save_rds must be either TRUE or FALSE")
        quit()
    }
}

#read in project file
# only take project data rows with a directory because some may be empty
project_data<-read.csv(project_data_path)
project_working<-project_data[project_data$directory!="",]

#get the path of the data produced by the post QC processing
seurat_object_path<-paste0(root,seurat_object_path)

if(!(
    (TRUE) &
    file.exists(seurat_object_path)
    ))
    {
    cat("driver check failed... One of the arguments is missing, or the RDS file specified does not exists.\n")
    cat(paste0("number of args: ",length(opt),"\n"))
    cat(paste0("seurat_object_exists: ",file.exists(seurat_object_path),"\n"))
    quit()
}

cat("driver check passed....\n")
source("annotate_cells_functions.R")

#create a folder to store results from manual annotation
subroot<-paste0(root,"annotate_cells_",dir_append,"/")
suppressWarnings(dir.create(subroot))
if (!(dir.exists(subroot))){
    cat("Results directory could not be created or could not be found. Quitting. \n")
    quit()
}else{
    cat("Results directory created \n")
}

#create the manual changes directory
suppressWarnings(dir.create(paste0(root,"manual_changes/")))

#read in the seurat object
cat("reading in the seurat file...\n")
seurat_object<-mcreadRDS(seurat_object_path,mc.cores=16)

if (use_manual_anno){
    #if we have run this analysis already maybe we want to modify clusters manually using a cluster modification file. Used to merge clusters together into one
    if(file.exists(paste0(root,"manual_changes/integrated_clusters_modification.tab"))){
        cluster_modify_info<-read_cluster_mod_data(paste0(root,"manual_changes/integrated_clusters_modification.tab"))
        seurat_object<-modify_clusters(seurat_object=seurat_object,cluster_mod_data=cluster_modify_info,subroot=subroot,root_external=root_external)
    }else{
        #if there is no integrated_clusters_modification.tab file then this is the first run, and we should create the file, so that the user can easily modify cluster afterwards in they want
        cluster_names<-data.frame(levels(seurat_object))
        colnames(cluster_names)<-"original"
        write.table(cluster_names,file=gdpath("integrated_clusters_modification.tab",root=paste0(root,"manual_changes/"),external=root_external),col.names=TRUE,row.names=FALSE,sep="\t")
    }
}else{
    cat("manual annotations not used.\n")
}

#only run deconvolution if specified to
if (run_deconvolution){
    #now lets add deconvolution info
    cat("\n##############################\n")
    cat("Adding scSNV deconvolution data\n")
    cat("##############################\n\n")
    all_recip_barcodes<-c()
    for (i in 1:nrow(project_working)){
        if (project_working$deconvolution_path[i]!="no_decon"){
            r_barcodes<-extract_recip_barcodes_from_gavin_table(path=project_working$deconvolution_path[i],label_in_label_column=project_working$decon_recip_label[i])
            cat(paste0("reading deconvolution data from: ",project_working$deconvolution_path[i],".",length(r_barcodes)," barcodes found in gavin's data.\n"))
            r_barcodes<-paste0(project_working$case[i],project_working$timepoint[i],"_",r_barcodes)
            all_recip_barcodes<-c(all_recip_barcodes,r_barcodes)
        }
    }
    #add metadata whether a cell is from donor or recipient
    seurat_object<-combined_add_deconvolution_info(seurat_object=seurat_object,recipient_barcodes=all_recip_barcodes)
    #also add the cell sex data for each cell 
    seurat_object<-add_cell_sex(seurat_object=seurat_object,project_data=project_data)    
    #doublet removal done in the background_decont combined script
}

#set idents to new ones based on manual anno
Idents(seurat_object)<-seurat_object@meta.data$integrated_clusters
print(levels(Idents(seurat_object)))
#set defulat assay to RNA for cell marker differential expression testing
DefaultAssay(seurat_object)<-"RNA"
seurat_object <- NormalizeData(seurat_object)


if(file.exists(paste0(root,"manual_changes/manual_cluster_names.tab"))){
    man_anno<-read.delim(paste0(root,"manual_changes/manual_cluster_names.tab"))
    if(length(levels(Idents(seurat_object)))==length(man_anno$man_anno)){
        cat("manual provided names are acceptable and will be used.\n")
        #full names
        manual_anno<-Idents(seurat_object)
        levels(manual_anno)<-man_anno$man_anno #later on Idents of cells get added as a metadata "GSVA clusters"
        seurat_object<-AddMetaData(object=seurat_object,metadata=manual_anno,col.name="GSVA_anno_labels")
        cat("added cell full names.\n")
        #nicknames
        manual_anno_nickname<-Idents(seurat_object)
        levels(manual_anno_nickname)<-man_anno$nickname
        seurat_object<-AddMetaData(object=seurat_object,metadata=manual_anno_nickname,col.name="GSVA_anno_labels_nickname")
        cat("added nicknames.\n")
    }else{
        cat("manual_cluster_names.tab file found but names are not acceptable.\n")
        cat(paste0("length of (levels(Idents(seurat_object)) is ",length(levels(Idents(seurat_object))),"\n"))
        cat(paste0("length of length(man_anno$man_anno) is ",length(man_anno$man_anno),"\n"))
        #annotate the cells using our premade gsva marker list
        annotate_cells<-annotate_combined_using_pre_only(
            seurat_object=seurat_object,
            infile_gmt=infile_gmt,
            root=subroot,
            root_external=root_external
        )
        manual_anno<-Idents(seurat_object)
        levels(manual_anno)<-annotate_cells
        seurat_object<-AddMetaData(object=seurat_object,metadata=manual_anno,col.name="GSVA_anno_labels")
    }
}else{
    #annotate the cells using our premade gsva marker list
    annotate_cells<-annotate_combined_using_pre_only(
        seurat_object=seurat_object,
        infile_gmt=infile_gmt,
        root=subroot,
        root_external=root_external
    )
    manual_anno<-Idents(seurat_object)
    levels(manual_anno)<-annotate_cells
    seurat_object<-AddMetaData(object=seurat_object,metadata=manual_anno,col.name="GSVA_anno_labels")
}

cat("\n############################################\n")
cat(paste0("there are ",length(unique(seurat_object$integrated_clusters))," clusters.\n"))
cat("############################################\n")

#create marker tables to annotate cells
if (make_marker_table){
    create_marker_table(
        seurat_object=seurat_object,
        min_prct_cells=min_prct_cells,
        timepoint="pre",
        logfc.threshold=logfc.threshold,
        root=subroot,
        root_external=root_external
    )

    create_marker_table(
        seurat_object=seurat_object,
        min_prct_cells=min_prct_cells,
        timepoint="post",
        logfc.threshold=logfc.threshold,
        root=subroot,
        root_external=root_external
    )
}

#create a pdf file with prelinary visualizations
pdf(gdpath("sc_non_integrated_integrated_clustering.pdf",root=subroot,external=root_external),height=10,width=17)
print(DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint","seurat_clusters","integrated_clusters","GSVA_anno_labels"),label=T,raster=FALSE))
print(DimPlot(seurat_object, reduction = "umap",group.by=c("integrated_clusters","GSVA_anno_labels"),label=T,raster=FALSE)+ NoLegend())
print(DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint"),label=T,raster=FALSE))
print(DimPlot(seurat_object, reduction = "umap",group.by=c("seurat_clusters"),label=T,raster=FALSE))
print(DimPlot(seurat_object, reduction = "umap",group.by=c("integrated_clusters"),label=T,raster=FALSE))
print(DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint","integrated_clusters"),label=T,raster=FALSE))
p2<-DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint"),label=F,raster=FALSE,combine=T) + ggtitle("Coloured by Timepoint") + NoLegend()
p1<-DimPlot(seurat_object, reduction = "umap",group.by=c("GSVA_anno_labels"),label=T,raster=FALSE,combine=T)+ ggtitle("Coloured by Cell Type")+ NoLegend()
grid.arrange(grobs=list(p1,p2),ncol=2)
print(DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint","seurat_clusters","integrated_clusters","GSVA_anno_labels"),label=F,raster=FALSE))
print(DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint"),label=T,raster=FALSE))
print(DimPlot(seurat_object, reduction = "umap",group.by=c("integrated_clusters"),label=T,raster=FALSE))
print(DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint","GSVA_anno_labels"),label=T,raster=FALSE))
#only make these figures if run_deconvolution true
if(run_deconvolution){
    print(DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint_recipient_origin"),label=T,raster=FALSE))
    print(DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint_recipient_origin","GSVA_anno_labels"),label=T,raster=FALSE))
    print(DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint_recipient_origin","GSVA_anno_labels"),label=F,raster=FALSE))
    p1<-DimPlot(seurat_object, reduction = "umap",group.by=c("GSVA_anno_labels"),label=T,raster=FALSE,combine=T)+ ggtitle("Coloured by Cell Type")+ NoLegend()
    p2<-DimPlot(seurat_object, reduction = "umap",group.by=c("timepoint_recipient_origin"),label=F,raster=FALSE,combine=T) + ggtitle("Coloured by Donor or Recipient") 
    grid.arrange(grobs=list(p1,p2),ncol=2)
    DimPlot(seurat_object, reduction = "umap",group.by=c("GSVA_anno_labels"),split.by="recipient_origin",label=F,raster=FALSE,combine=T) + ggtitle("Coloured by Donor or Recipient")
}
print(DimPlot(seurat_object, reduction = "umap",cells=seurat_object$timepoint=="pre",group.by=c("GSVA_anno_labels"),label=T,raster=FALSE))
print(DimPlot(seurat_object, reduction = "umap",cells=seurat_object$timepoint=="post",group.by=c("GSVA_anno_labels"),label=T,raster=FALSE))
print(DimPlot(seurat_object, reduction = "umap",group.by=c("GSVA_anno_labels"),label=T,raster=FALSE))
for (case in unique(project_working$case)){
    print(DimPlot(seurat_object, reduction = "umap",cells=seurat_object@meta.data$ltx_case==case,group.by=c("timepoint"),label=T,raster=FALSE) + plot_annotation(title = paste0(case," by timepoints.")))
}
dev.off()


if(make_save_rds){
        #save a copy which contains all cells
        cat("Saving a copy of the seurat object with recipient cells.\n")
        mcsaveRDS(seurat_object,file=gdpath(paste0(project_name,"_man_anno_WITH_recipient_cells.rds"),root=subroot,external=root_external),mc.cores=16)
    }

#lets remove recipient cells, do that only donor cells remain.
if(run_deconvolution){
    cat("saving seurat object with only donor cells\n")
    seurat_object<-seurat_object[,seurat_object@meta.data$recipient_origin==FALSE]
}

#save output
if(make_save_rds){
    mcsaveRDS(seurat_object,file=gdpath(paste0(project_name,"_man_anno.rds"),root=subroot,external=root_external),mc.cores=16)
}
cat("complete.\n")

#save a text summary of cells per case
sink(gdpath("integratedClusters_ltxCase_breakdown.txt",root=subroot,external=root_external))
table(seurat_object$integrated_clusters,seurat_object$ltx_case)
sink()