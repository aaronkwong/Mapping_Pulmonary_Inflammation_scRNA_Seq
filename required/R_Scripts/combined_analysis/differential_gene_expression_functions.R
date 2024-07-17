library(trqwe)
library(Seurat)
library(doParallel)
registerDoParallel(cores=6)
library(limma)
library(pcaPP)
library(scales)


calc_mean_logExpression<-function(x,pseudocount.use=1){
	return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
}

calc_log2FC<-function(x,pseudocount.use=1){
	return(log2(x = mean(x = expm1(x = x)) + pseudocount.use))
}


library(nebula)
#nebula incorporates a random effect model to account for bias of pseudoreplication when doing DE
run_nebula_raw_counts<-function(seurat_object,metadata,cell_cluster,cpc=0.05,min_unique_genes_per_cell=200){
	cat(paste0("nebula: processing cell type ",cell_cluster,"\n"))
	x0<-seurat_object[,seurat_object$integrated_clusters==cell_cluster]
	#lets get the rna raw counts
	x0_counts<-x0@assays$RNA@counts
	x0_counts<-x0_counts[!apply(x0_counts,1,FUN=function(x){all(x==0)}),] #remove any genes where all expression in 0
	# cells_removed<-!apply(x0_counts,2,FUN=function(x){all(x==0)}) #remove any cells where library size is 0
	cells_removed<-apply(x0_counts,2,FUN=function(x){sum(x!=0)})>min_unique_genes_per_cell # for selecting cells with more than 200 unique genes
	x0_counts<-x0_counts[,cells_removed]
	#lets get sample IDs
	x0_sid<-paste0(x0@meta.data$timepoint,x0@meta.data$ltx_case)
	x0_sid<-x0_sid[cells_removed]
	#create design matrix
	pred<-x0@meta.data[,c("timepoint","ltx_case")]
	pred<-pred[cells_removed,]
	pred$timepoint[pred$timepoint=="pre"]<-0
	pred$timepoint[pred$timepoint=="post"]<-1
	df = model.matrix(~timepoint+ltx_case, data=pred)
	#norm factors
	norm_factors<-colSums(as.matrix(x0_counts))
	# norm_factors<-calcNormFactors(x0_counts,method = c("TMM"))
	re = nebula(x0_counts,x0_sid,pred=df,cpc=cpc,offset=norm_factors)
	summ<-re$summary
	summ$convergence<-re$convergence
	summ<-summ[order(summ$p_timepoint1,decreasing=FALSE),]
	summ$p_val_adj<-p.adjust(summ$p_timepoint1,method="fdr")
	summ_final<-summ[,c("gene","logFC_timepoint1","p_timepoint1","p_val_adj","convergence")]
	rownames(summ_final)<-summ_final$gene
	colnames(summ_final)<-c("gene","avg_logFC","p_val","p_val_adj","convergence")
	#count pct expression 
	ident.1.pct<-apply(x0_counts,MARGIN=1,FUN=function(expr,group){return((sum(expr[group==1]!=0)/sum(group==1)))},group=pred$timepoint) #pred$recipient_origin=0, 0 is pre
	#count number of cells expressing at least 1 copy of this transcript in group 2
	ident.2.pct<-apply(x0_counts,MARGIN=1,FUN=function(expr,group){return((sum(expr[group==0]!=0)/sum(group==0)))},group=pred$timepoint) #pred$recipient_origin=1, 1 is post
	#add these to the final table
	summ_final$ident.1.pct<-ident.1.pct[rownames(summ_final)]
	summ_final$ident.2.pct<-ident.2.pct[rownames(summ_final)]
	return(summ_final)
}

#this function does differential gene expression between ident.1 and ident.2, intended to compare donor vs recipient
# ident.1.pct is donor pct
# ident.2.pct is recip pct
run_nebula_raw_counts_deconvolution<-function(seurat_object,ident.1, ident.2,metadata,cpc=0.05,min_unique_genes_per_cell=200){
	cat(paste0("nebula: comparing ident.1 to ident.2 ","\n"))
	x0<-seurat_object[,c(ident.1, ident.2)]
	print(table(x0$recipient_origin,x0$ltx_case))
	#lets get the raw counts
	x0_counts<-x0@assays$RNA@counts
	x0_counts<-x0_counts[!apply(x0_counts,1,FUN=function(x){all(x==0)}),] #remove any genes where all expression is 0
	# cells_removed<-!apply(x0_counts,2,FUN=function(x){all(x==0)}) #remove any cells where library size is 0
	cells_removed<-apply(x0_counts,2,FUN=function(x){sum(x!=0)})>min_unique_genes_per_cell # for selecting cells with more than 200 unique genes
	x0_counts<-x0_counts[,cells_removed]
	#lets get sample IDs
	x0_sid<-paste0(x0@meta.data$ltx_case,x0@meta.data$recipient_origin)
	x0_sid<-x0_sid[cells_removed]
	#create design matrix
	pred<-x0@meta.data[,c("recipient_origin","ltx_case","cell_sex")]
	pred$recipient_origin<-as.character(pred$recipient_origin)
	pred$recipient_origin[pred$recipient_origin=="TRUE"]<-0 #this will make donor/recip comparison
	pred$recipient_origin[pred$recipient_origin=="FALSE"]<-1
	pred<-pred[cells_removed,]
	df = model.matrix(~recipient_origin+ltx_case, data=pred)
	#norm factors
	norm_factors<-colSums(as.matrix(x0_counts))
	# norm_factors<-calcNormFactors(x0_counts,method = c("TMM"))
	re = nebula(x0_counts,x0_sid,pred=df,cpc=cpc,offset=norm_factors)
	summ<-re$summary
	summ$convergence<-re$convergence
	summ<-summ[order(summ$p_recipient_origin1,decreasing=FALSE),]
	summ$p_val_adj<-p.adjust(summ$p_recipient_origin1,method="fdr")
	summ_final<-summ[,c("gene","logFC_recipient_origin1","p_recipient_origin1","p_val_adj","convergence")]
	rownames(summ_final)<-summ_final$gene
	colnames(summ_final)<-c("gene","avg_logFC","p_val","p_val_adj","convergence")
	#count pct expression 
	ident.1.pct<-apply(x0_counts,MARGIN=1,FUN=function(expr,group){return((sum(expr[group==1]!=0)/sum(group==1)))},group=pred$recipient_origin) #pred$recipient_origin=0, 0 is the donor
	#count number of cells expressing at least 1 copy of this transcript in group 2
	ident.2.pct<-apply(x0_counts,MARGIN=1,FUN=function(expr,group){return((sum(expr[group==0]!=0)/sum(group==0)))},group=pred$recipient_origin) #pred$recipient_origin=1, 1 is the recip
	#add these to the final table
	summ_final$ident.1.pct<-ident.1.pct[rownames(summ_final)]
	summ_final$ident.2.pct<-ident.2.pct[rownames(summ_final)]
	return(summ_final)
}