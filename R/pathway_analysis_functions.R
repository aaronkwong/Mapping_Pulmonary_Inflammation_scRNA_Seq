library(trqwe)
library(Seurat)
library(doParallel)
registerDoParallel(cores=6)
library(limma)
library(pcaPP)
library(scales)
library(nebula)



# this function computes differential gene expression between ident.2 vs ident.1
# assumed that seurat objects contains only cells that are part of comparison. Should throw an error if a cell not being compared in included
# use case data as fixed effects
# used for deconvolution pathway analysis
# can handle single cell data containg cells from multiple cases or just a single case
run_BPSC_de<-function(seurat_object,ident.1,ident.2,pseudocount.use=1){
	cat(paste0("BPSC: extracting count matrix\n"))
	# if (!all(rownames(seurat_object@meta.data) %in% append(ident.1,ident.2))){
	# 	cat("Not all barcodes in seurat object found in append(ident.1,ident.2). Make sure that all barcodes provided in ident.1 and ident.2 cover all barcodes found in the seurat object.\n")
	# }
	clus<-seurat_object[,rownames(seurat_object@meta.data) %in% append(ident.1,ident.2)]
	#extract counts
	clus_expr<-clus@assays$SCT@counts
	#remove any rows that are all NAs
	non_na_index<-apply(clus_expr,MARGIN=1,FUN=function(x){!all(is.na(x))})
	clus_expr<-clus_expr[non_na_index,]
	#get vector for blocking representing timepoint
	group<-rownames(clus@meta.data)
	group[group %in% ident.1]<-1
	group[group %in% ident.2]<-2
	group<-as.numeric(group)
	#get vector for blocking represenitng ltx_case
	case<-clus@meta.data$ltx_case
	#set controls 
	controlIds=which(group==1)
	#Create a design matrix including the group labels. All batch effects can be also added here if they are available
	#check if multiple cases are present in our cells for analysis, and set cases as fixed effects if needed.
	if(length(unique(case))>1){
		design=model.matrix(~group+case) #if we are grouping multiple cases together we can set each case as a fixed effect
	}else{
		design=model.matrix(~group) #if we are doing per-case analysis then we dont need case fixed effects
	}
	cat(paste0("BPSC: fitting model\n"))
	#Select the column in the design matrix corresponding to the coefficient (the group label) for the GLM model testing
	coef=2 
	try_class<-try(res<-BPglm(data=RelativeCounts(clus_expr, scale.factor = 1e6, verbose = TRUE), controlIds=controlIds, design=design, coef=coef, estIntPar=FALSE, useParallel=TRUE))
	if(class(try_class) %in% "try-error"){
		return("error")
	}
	gc()
	cat(paste0("BPSC: calculating FC and pct expression\n"))
	#now we need to compute fold change as BPSC doesnt compute. Lets use logFC
	# we need to pull log transformed data
	clus_expr_log<-clus@assays$SCT@data
	#REMOVE ANY ROWS using the same index used to filter count data (or else there will be different number of rows between BPglm and logFC results)
	clus_expr_log<-clus_expr_log[non_na_index,]
	#we need to check to make sure these logFC are accurate
	res$logFC<-apply(clus_expr_log,MARGIN=1,FUN=function(log_expr,group){return(calc_mean_logExpression(log_expr[group==2])-calc_mean_logExpression(log_expr[group==1]))},group=group)
	#count number of cells expressing at least 1 copy of this transcript in group 1
	res$ident.1.pct<-apply(clus_expr_log,MARGIN=1,FUN=function(log_expr,group){return((sum(log_expr[group==1]!=0)/sum(group==1)))},group=group)
	#count number of cells expressing at least 1 copy of this transcript in group 2
	res$ident.2.pct<-apply(clus_expr_log,MARGIN=1,FUN=function(log_expr,group){return((sum(log_expr[group==2]!=0)/sum(group==2)))},group=group)
	res$gene_scores<--log(res$PVAL)*sign(res$logFC)
	res$ranked_list_desc<-names(res$PVAL)[order(res$gene_scores,decreasing=TRUE)]
	return(res)
}

#this function is a wrapper for limma_voom. 
# this function computes differential gene expression between ident.2 vs ident.1
# assumed that seurat objects contains only cells that are part of comparison. Should throw an error if a cell not being compared in included
# use case data as fixed effects
# used for deconvolution pathway analysis
# can handle single cell data containg cells from multiple cases or just a single case
run_de_limma_voom<-function(seurat_object,ident.1,ident.2){
	#seprate object into the cell type we want
	if(!(append(ident.1,ident.2) %in% rownames(seurat_object@meta.data))){
		cat("at least one ident (1 or 2) is not present in the suerat object")
	}
	clus<-seurat_object[,rownames(seurat_object@meta.data) %in% append(ident.1,ident.2)]
	#pull raw expression matrix
	clus_expr<-clus@assays$RNA@counts
	case<-clus@meta.data$ltx_case
	if(length(unique(case))>1){
		return(limma_voom(count_matrix=clus_expr,ident.1=ident.1,ident.2=ident.2,case=case))
	}else{
		return(limma_voom(count_matrix=clus_expr,ident.1=ident.1,ident.2=ident.2))
	}

}


# built using example https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_voomlimma.R
# https://support.bioconductor.org/p/85511/
# https://f1000research.com/articles/5-1408
# this is the basic limma voom function used by run_cell_type_de_limma_voom
limma_voom<-function(count_matrix,ident.1,ident.2,case=NULL,pct.threshold=0.1){
	#remove any rows where expression is all zero
	counts<-count_matrix[!apply(count_matrix,1,FUN=function(x){all(x==0)}),]

	if(sum(sum(colnames(counts) %in% ident.1),sum(colnames(counts) %in% ident.2)) == sum(length(ident.1),length(ident.2)))
	{
		group<-rep(0,length(ident.1)+length(ident.2))
		group[colnames(counts) %in% ident.1]<-0
		group[colnames(counts) %in% ident.2]<-1
	}else{
		cat("error, the idents provided do not perfectly match colnames in the count matrix")
	}

	#get vector for blocking representing ltx_case
	if (!is.null(case)){
		case<-factor(case)
		design=model.matrix(~group+case)
	}else{
		design=model.matrix(~group)
	}
	print(design)

	#filtering out lowly expressed genes improves voom-limma performance. Here we will filter based on pct cells expressing like in seurat
	#count pct expression 
	ident.1.pct<-apply(counts,MARGIN=1,FUN=function(expr,group){return((sum(expr[group==0]!=0)/sum(group==0)))},group=group)
	#count number of cells expressing at least 1 copy of this transcript in group 2
	ident.2.pct<-apply(counts,MARGIN=1,FUN=function(expr,group){return((sum(expr[group==1]!=0)/sum(group==1)))},group=group)
	#remove any rows of lowly expressed genes
	counts<-counts[(ident.1.pct>pct.threshold | ident.2.pct>pct.threshold),]
	ident.1<-ident.1.pct[(ident.1.pct>pct.threshold | ident.2.pct>pct.threshold)]
	ident.2<-ident.2.pct[(ident.1.pct>pct.threshold | ident.2.pct>pct.threshold)]
	cat(paste0(length(ident.1)," genes passed expression threshold for differential gene expression analysis.\n"))

	cat("limma voom: calculating normalization factors\n")
	dge<-DGEList(counts)
	dge<- calcNormFactors(dge)
	
	# voom
	cat("limma voom: fitting model\n")
	vm <- voom(dge, design = design, plot = TRUE)
	fit <- lmFit(vm, design = design)
	fit <- eBayes(fit)
	tt <- topTable(fit,coef="group", n = Inf, adjust.method = "BH")
	#count pct expression 
	tt$ident.1.pct<-ident.1
	#count number of cells expressing at least 1 copy of this transcript in group 2
	tt$ident.2.pct<-ident.2
	colnames(tt)[colnames(tt)=="P.Value"]<-"p_val"
	colnames(tt)[colnames(tt)=="logFC"]<-"avg_logFC"
	colnames(tt)[colnames(tt)=="adj.P.Val"]<-"p_val_adj"
	return(tt)
}


#run cell type differential gene expression
run_cell_type_de_seurat_LR<-function(seurat_object,metadata,cell_cluster,thresh=0.1,logFC_thresh=0){
	cat(paste0("Seurat_LR: analyzing cell type ",cell_cluster,"\n"))
	#seprate object into the cell type we want
	clus<-seurat_object[,seurat_object@meta.data[,metadata]==cell_cluster]
	if(length(unique(clus$ltx_case))>1){
		t<-FindMarkers(clus,ident.1=rownames(clus@meta.data)[clus$timepoint=="post"],ident.2=rownames(clus@meta.data)[clus$timepoint=="pre"],test.use="LR",latent.vars = "ltx_case",min.pct = thresh,logfc.threshold=logFC_thresh)
	}else{
		t<-FindMarkers(clus,ident.1=rownames(clus@meta.data)[clus$timepoint=="post"],ident.2=rownames(clus@meta.data)[clus$timepoint=="pre"],test.use="LR",min.pct = thresh,logfc.threshold=logFC_thresh)
	}
	colnames(t)[colnames(t)=="avg_log2FC"]<-"avg_logFC"
	return(t)
}

#negbinom
run_cell_type_de_seurat_negbinom<-function(seurat_object,metadata,cell_cluster,thresh=0.1,logFC_thresh=0){
	cat(paste0("Seurat_negbinom: analyzing cell type ",cell_cluster,"\n"))
	#seprate object into the cell type we want
	clus<-seurat_object[,seurat_object@meta.data[,metadata]==cell_cluster]
	if(length(unique(clus$ltx_case))>1){
		t<-FindMarkers(clus,ident.1=rownames(clus@meta.data)[clus$timepoint=="post"],ident.2=rownames(clus@meta.data)[clus$timepoint=="pre"],test.use="negbinom",latent.vars = "ltx_case",min.pct = thresh,logfc.threshold=logFC_thresh)
	}else{
		t<-FindMarkers(clus,ident.1=rownames(clus@meta.data)[clus$timepoint=="post"],ident.2=rownames(clus@meta.data)[clus$timepoint=="pre"],test.use="negbinom",min.pct = thresh,logfc.threshold=logFC_thresh)
	}
	colnames(t)[colnames(t)=="avg_log2FC"]<-"avg_logFC"
	return(t)
}

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

# This function takes a list, where each object in the list is a list of differenitally expresssed genes per cluster and conducts enrichement analysis.
# @param rank_files list, where each element is a list of DEG
# @param min_geneset geneset min size
# @param max_geneset geneset max size
# @param gmt_database gmt database
# @param out_folder output directory
# @param root The root directory
# @param root_external whether the root is in the cloud synced folder
# @param project_name The name of the project
# @param output_append prefix to append to results directory

run_GSEA<-function(rank_files,min_geneset,max_geneset,gmt_database,out_folder,root,root_external,project_name,output_append="enrichment"){
	#here we will correlate each item name "cluster_1" with the path "path/to/cluster_1.rnk"
	for (i in 1:(length(rank_files))){
		if(rank_files[i]!="insufficient_cells"){
			sent<-0
			while(sent!=1){
				#you might chnage this value if the computer does not have much ram
				if (available_ram()>50000){
					s.broad_GSEA_preranked(
					s.gsea_soft=gdpath("Masters/Software/gsea-3.0.jar"),
					s.gsea_memory=4096,
					s.gsea_nperm=1000,
					s.gsea_rnk=short_path(rank_files[i]),
					s.gsea_gmt=short_path(gmt_database),
					s.gsea_output_name=paste0(names(rank_files)[i],"_",output_append,"_"),
					s.gsea_output_location=out_folder, #this cannot be in Onedrieve due to hyphens in the name
					s.timestamp=123456789,
					s.set_min=min_geneset,
					s.set_max=max_geneset,
					wait=FALSE #only do this if you have the ram to run all enrichment concurrently
					)
					print(paste0("analysis for cluster ",i," has been sent."))
					Sys.sleep(15)
					sent<-1
				}else{
					cat("insufficient memory available... please free up ram\n")
					Sys.sleep(1)
				}
			}
		}else{
			cat("no DE to run due to insufficient cells.\n")
		}
	}
}


# @param seurat_object seurat object containing cells
# @param subroot directory for results
# @param root The root directory
# @param root_external whether the root is in the cloud synced folder
# @param thresh is a number between 0 and 1 defining the minimum proportion of cells that must be expressed in either ident.1 or ident.2 to be kept for analysis
# @param logFC_thresh logFC threshold for genes to keep
# @param de_tool the method for differential gene expression testing. Can be nebula, voomlimma, BPSC, seurat_neg_binom, seurat_LR
run_pathway_analysis_combined<-function(seurat_object,subroot,root,root_external,project_name,thresh=0,logFC_thresh=0,de_tool="nebula"){
	if(file.exists(gdpath(paste0(project_name,"_all_cluster_DE.rds"),root=subroot,external=root_external))){
		df.cluster<-readRDS(file=gdpath(paste0(project_name,"_all_cluster_DE.rds"),root=subroot,external=root_external))
		rank_files<-readRDS(file=gdpath(paste0(project_name,"DE_rank_files.rds"),root=subroot,external=root_external))
	}else{

		cat("creating pathway_analysis directory output...\n")
		dir.create(subroot)
		cat("starting pathway analysis...\n")
		#create the compound label so that cells of the same cluster but are from different timepoints are labelled accordingly
		compound_labels<-paste0(seurat_object$timepoint,seurat_object$integrated_clusters)
		seurat_object<-AddMetaData(object=seurat_object,metadata=compound_labels,col.name="compound")
		Idents(seurat_object)<-seurat_object$compound
		#lets do a check to make sure that there are enough cells pre and post in each cluster
		if (all(table(seurat_object$timepoint,seurat_object$integrated_clusters)<5)){
			cat("\nError, some clusters have fewer than 5 cells at at least one timepoint: \n")
			print(table(seurat_object$timepoint,seurat_object$integrated_clusters))
			quit()
		}
		#compute differential gene expression for each cluster
		#create vector to store rank file locations
		cat("Calculating cluster specific differential gene expression and making p-value histogram plots...\n")
		pdf(gdpath(paste0("cluster_pathway_analysis_",project_name,".pdf"),root=subroot,external=root_external))
		rank_files<-c()
		df.cluster<-list()
		# df.cluster_name<-c()
		cat(paste0("There are ",length(levels(seurat_object@meta.data$integrated_clusters))," cell types detected.\n"))
		for (i in 0:(length(levels(seurat_object@meta.data$integrated_clusters))-1)) {
		# for (i in 0) {
			#make sure there are enough pre and post cells
			num_pre<-sum(seurat_object@meta.data$compound==paste0("pre",i))
			num_post<-sum(seurat_object@meta.data$compound==paste0("post",i))
			if(num_pre>5 & num_post>5){#if there are enough cells for DE
				if(de_tool=="BPSC"){
					#logfc limit does not benefit by change to -Inf due 
					test.df <- run_cell_type_de(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i,pseudocount.use=1)
					
					#1 must be added since counter starts from 0
					
					de_res<-data.frame(test.df$PVAL,test.df$logFC,p.adjust(test.df$PVAL,method="BH"),test.df$ident.1.pct,test.df$ident.2.pct)
					colnames(de_res)<-c("","avg_logFC","p_val_adj","ident.1.pct","ident.2.pct")
					rownames(de_res)<-names(test.df$PVAL)
					de_res<-de_res[!is.na(de_res$p_val),]
					de_res<-de_res[order(de_res$p_val_adj),]
					de_res<-de_res[(de_res$ident.1.pct>thresh | de_res$ident.2.pct>thresh),]
				}else if(de_tool=="voomlimma"){ #limma voom method
					de_res<-run_cell_type_de_limma_voom(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i)
				}else if(de_tool=="seurat_LR"){ #seurat LR
					de_res<-run_cell_type_de_seurat_LR(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i,thresh=thresh,logFC_thresh=logFC_thresh)
				}else if(de_tool=="seurat_negbinom"){ #seurat negbinom
					de_res<-run_cell_type_de_seurat_negbinom(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i,thresh=thresh,logFC_thresh=logFC_thresh)
				}else if(de_tool=="nebula"){ #seurat negbinom
					de_res<-run_nebula_raw_counts(seurat_object=seurat_object,metadata="integrated_clusters",cell_cluster=i)
				}
				df.cluster[[i+1]]<-de_res

				#make historgram
				print(hist(de_res$p_val,breaks=20,main=paste0("cluster_",i,"_min_express_cutoff_",thresh)))

				#make data frame for make_rnk
				# df<-data.frame(de_res$p_val,de_res$avg_logFC)
				# df<-df[!is.na(df$p_val) & !is.na(df$avg_logFC),]
				
				#FIX P VALUES SO SMALL THEY ARE REPORTED AS -Inf
				#some p values are so small that they are reported as 0. We need to change these values to something that reflects their position for the ranked list. 
				#make_rnk_simple(test.df[,1:2],paste0("C:/Users/Aaron Wong/OneDrive - UHN/Masters/project_6_scRNA_seq/analyze_sc/practice/de_gsea/cluster_",i,".rnk"),clean=F)
				# cat("using make_rnk_simple_by_fold_change\n")
				# make_rnk_simple_by_fold_change(de_res[,c("p_val","avg_logFC")],gdpath(paste0(project_name,"_cluster_",i,".rnk"),root=subroot,external=root_external),clean=F)
				# make_rnk_simple(de_res[,c("p_val","avg_logFC")],gdpath(paste0(project_name,"_cluster_",i,".rnk"),root=subroot,external=root_external),clean=F)
				rank_files<-c(rank_files,gdpath(paste0(project_name,"_cluster_",i,".rnk"),root=subroot,external=root_external))
			}else{
				df.cluster[[i+1]]<-"insufficient_cells"
				rank_files<-c(rank_files,"insufficient_cells")
			}

			# df.cluster_name<-c(df.cluster_name,paste0("cluster_",i))
			print(paste0(i," run complete"))
		}
		dev.off()
		names(rank_files)<-paste0(project_name,"_cluster_",0:as.numeric(length(rank_files)-1))
		names(df.cluster)<-paste0("cluster_",0:as.numeric(length(rank_files)-1))
		#df.cluster is a list object where each list corresponds to the differnetially expressed genes for each cluster
		#rank.files is a list of the names of the rank files generated
		saveRDS(df.cluster,file=gdpath(paste0(project_name,"_all_cluster_DE.rds"),root=subroot,external=root_external))
		saveRDS(rank_files,file=gdpath(paste0(project_name,"DE_rank_files.rds"),root=subroot,external=root_external))

		#save cluster DE
		cat("Saving cluster differential expression....\n")
		for (i in 1:length(df.cluster)){
			write.table(df.cluster[[i]],file=gdpath(paste0(project_name,"_cluster_",i-1,"_DE.tab"),root=subroot,external=root_external),sep="\t",row.names=TRUE,col.names=NA)
		}
		cat("Pathway analysis complete....\n")
	}

	#cycle through cell numbers and grab the DE tables. Then cutoff by thresh (thresh is the minimum proportion of cells in either ident.1 or ident2 that must be met for a gene to be kept)
	for (i in 0:(length(levels(seurat_object@meta.data$integrated_clusters))-1)) {
		if(df.cluster[[i+1]]!="insufficient_cells"){
			subset<-df.cluster[[i+1]][(df.cluster[[i+1]][,"ident.1.pct"]>thresh | df.cluster[[i+1]][,"ident.2.pct"]>thresh),]
			cat("using make_rnk_simple_by_fold_change\n")
			#current results based on this Nov 3, 2021
			# make_rnk_simple_by_fold_change(subset[,c("p_val","avg_logFC")],gdpath(paste0(project_name,"_cluster_",i,".rnk"),root=subroot,external=root_external),clean=F)
			make_rnk_simple(subset[,c("p_val","avg_logFC")],gdpath(paste0(project_name,"_cluster_",i,".rnk"),root=subroot,external=root_external),clean=F)
		}else{
			cat(paste0("cluster ",i," had insufficient cells.\n"))
		}
	}

	return(rank_files)
}


