# This script takes a seurat object containing recipient cells and removes recipient cells from clusters where they are not expected. Generates plots and summary tables of recipient cells.

suppressPackageStartupMessages(library(renv))
renv::activate()
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(trqwe))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
source("gdpath.R")

option_list <- list(
    make_option(c("-p", "--seurat_object_path"), default="NA",
        help="The seurat object"),
    #
    make_option(c("-e", "--clusters_with_recip_cells"), default="none",
        help="cluster numbers containing recip cells seperated by commas"),
    #
    make_option(c("-v", "--root"), default="NA",
        help="The directory where the various result directories for each part of the pipeline will be written to."),
    #
    make_option(c("-q", "--root_external"), default="NA",
        help="TRUE or FALSE. If the root is contatined in the cloud synced folder."),
    #
    make_option(c("-a", "--make_save"), default=FALSE,
        help="")
  )

opt <- parse_args(OptionParser(option_list=option_list))

#set params
seurat_object_path<-opt$seurat_object_path
root<-opt$root
if(opt$make_save=="TRUE"){
    make_save<-TRUE
}else{
    make_save<-FALSE
}

#set params for output
if(opt$root_external=="TRUE"){
    root_external<-TRUE
}else{
    root_external<-FALSE
    cat("internal root detected.\n")
}
raw_recip_data<-opt$clusters_with_recip_cells

#create output directory
subroot<-paste0(root,"deconvolution_analysis/")
suppressWarnings(dir.create(subroot))

# load seurat object containing cells (and the recipient cells)
seurat_object<-mcreadRDS(seurat_object_path)

#plot all donor and recip cells form all cases
pdf(gdpath("deconvolution_visualization.pdf",root=subroot,external=root_external),height=10,width=16)
    DimPlot(seurat_object,group.by=c("recipient_origin"),order=c("TRUE","FALSE"),raster=FALSE)+ plot_annotation(title = paste0("all cases donor and recip breakdown"))
    DimPlot(seurat_object,group.by=c("recipient_origin"),order=c("TRUE","FALSE"),raster=FALSE, label=TRUE)+ plot_annotation(title = paste0("all cases donor and recip breakdown"))
    #plot breakdown of recip and donor cells per case
    for (case in unique(seurat_object$ltx_case)){
    	print(DimPlot(seurat_object,cells=seurat_object$ltx_case==case,group.by=c("recipient_origin","integrated_clusters"))+ plot_annotation(title = paste0(case," donor and recip breakdown")))
    }
dev.off()

if (length(raw_recip_data)>1){
    run_recip_only_vis<-TRUE
}else{
    if(raw_recip_data!="none"){
        run_recip_only_vis<-TRUE
    }else{
        run_recip_only_vis<-FALSE
    }
}

#remove the small number random recipient parenchymal cells, by only keeping recipient immune cell clusters specified
if(run_recip_only_vis){
    recipient_clusters<-as.numeric(unlist(strsplit(raw_recip_data,split=",")))
    seurat_object_only_expected_recipient_cells<-seurat_object[,!(!(seurat_object$integrated_clusters %in% recipient_clusters) & seurat_object$recipient_origin==TRUE)]
    pdf(gdpath("deconvolution_visualization_only_expected_recip_cells.pdf",root=subroot,external=root_external),height=10,width=16)
        p1<-DimPlot(seurat_object_only_expected_recipient_cells, reduction = "umap",group.by=c("GSVA_anno_labels"),label=T,raster=FALSE,combine=T)+ ggtitle("Coloured by Cell Type")+ NoLegend()
        p2<-DimPlot(seurat_object_only_expected_recipient_cells, reduction = "umap",group.by=c("timepoint_recipient_origin"),label=F,combine=T,raster=FALSE) + ggtitle("Coloured by Donor or Recipient")
        grid.arrange(grobs=list(p1,p2),ncol=2)
        print(DimPlot(seurat_object_only_expected_recipient_cells,split.by="recipient_origin",group.by="GSVA_anno_labels",raster=FALSE,label=TRUE))
    dev.off()
}

# write out the seurat object after removing recipient cells that clustered with parenchymal cells
if(make_save){
  mcsaveRDS(seurat_object_only_expected_recipient_cells,file=gdpath(paste0("seurat_object_only_expected_recipient_cells",".rds"),root=subroot,external=root_external))
}

#output visualization of cells long name version
png(gdpath("deconvolution_visualization_only_expected_recip_cells.png",root=subroot,external=root_external),width=1200,height=900)
print(DimPlot(seurat_object_only_expected_recipient_cells,split.by="recipient_origin",group.by="GSVA_anno_labels",raster=FALSE,label=TRUE,label.size=8))
dev.off()

#output visualization of cells nickname version
png(gdpath("deconvolution_visualization_only_expected_recip_cells_nickname.png",root=subroot,external=root_external),width=1200,height=900)
print(DimPlot(seurat_object_only_expected_recipient_cells,split.by="recipient_origin",group.by="GSVA_anno_labels_nickname",raster=FALSE,label=TRUE,label.size=8,repel=T))
dev.off()

#output summary with counts of all cells, by cell type per case
sink(gdpath("WHOLE_DATASET_cell_type_breakdown_by_case_report.txt",root=subroot,external=root_external))
  table(seurat_object$integrated_clusters,seurat_object$ltx_case)
sink()

#output summary with counts of cells from the POST transplant timepoint, by cell type per case
seurat_object_post<-seurat_object[,seurat_object$timepoint=="post"]
sink(gdpath("POST_ONLY_donor_recipient_breakdown_report.txt",root=subroot,external=root_external))
	table(seurat_object_post$integrated_clusters,seurat_object_post$recipient_origin,seurat_object_post$ltx_case)
sink()

#same as one above except that the 3d table converted to a 2 dimensional table
sink(gdpath("POST_ONLY_donor_recipient_breakdown_report_flattened.txt",root=subroot,external=root_external))
    a<-table(seurat_object_post$integrated_clusters,seurat_object_post$recipient_origin,seurat_object_post$ltx_case)
    cases<-attributes(a)$dimnames[[3]]
    case_colnames<-c()
    for (case in cases){
        case_colnames<-c(case_colnames,paste0(rep(case,2),c("_Recip_FALSE","_Recip_TRUE")))
    }
    dim(a)<-c(dim(a)[1],dim(a)[2]*dim(a)[3])
    rownames(a)<-seq(0,nrow(a)-1,1)
    colnames(a)<-case_colnames
sink()

#output summary with counts of cells from the PRE transplant timepoint, by cell type per case
seurat_object_pre<-seurat_object[,seurat_object$timepoint=="pre"]
sink(gdpath("PRE_ONLY_breakdown_report.txt",root=subroot,external=root_external))
    table(seurat_object_pre$integrated_clusters,seurat_object_pre$ltx_case)
sink()