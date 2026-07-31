#Sections that write files are comented
#Loading packages
#Some packages belong to Bioconductor, to install them, see https://www.bioconductor.org/
#Package versions for all three scripts: doRNG_1.8.6.3, rngtools_1.5.2, doParallel_1.0.17, iterators_1.0.14, foreach_1.5.2, circlize_0.4.18, ComplexHeatmap_2.24.1, msigdbr_25.1.1, writexl_1.5.4, limma_3.64.3, ActivePathways_2.0.6  GSEMA_0.99.4, biomaRt_2.64.0, mathjaxr_2.0-0,     RColorBrewer_1.1-3, rstudioapi_0.18.0,  jsonlite_2.0.0,     shape_1.4.6.1,      magrittr_2.0.4,     magick_2.9.0,       farver_2.1.2,       GlobalOptions_0.1.4, vctrs_0.7.1,, memoise_2.0.1,      S4Arrays_1.8.1,     progress_1.2.3,     curl_7.0.0,, Rhdf5lib_1.30.0,    SparseArray_1.8.1,  rhdf5_2.52.1, plyr_1.8.9,, httr2_1.2.2,, impute_1.82.0,      cachem_1.1.0,       lifecycle_1.0.5,    pkgconfig_2.0.3,    rsvd_1.0.5, Matrix_1.7-4,       R6_2.6.1, , fastmap_1.2.0,      GenomeInfoDbData_1.2.14     rbibutils_2.4.1,    MatrixGenerics_1.20.0       clue_0.3-67,, digest_0.6.39,      numDeriv_2016.8-1.1, colorspace_2.1-2, AnnotationDbi_1.70.0, S4Vectors_0.46.0,   irlba_2.3.7,, pkgload_1.5.0,      GenomicRanges_1.60.0, RSQLite_2.4.6,      beachmat_2.24.0, filelock_1.0.3,     metadat_1.4-0,      httr_1.4.8,, abind_1.4-8,, compiler_4.5.2,     bit64_4.6.0-1,      S7_0.2.1, , metafor_4.8-0,      BiocParallel_1.42.1, DBI_1.2.3, HDF5Array_1.36.0,   rappdirs_0.3.4,     DelayedArray_0.34.1, rjson_0.2.23,       tools_4.5.2,, otel_0.2.0,, glue_1.8.0,, h5mread_1.0.1,      nlme_3.1-168,       rhdf5filters_1.20.0, cluster_2.1.8.1,    generics_0.1.4,     gtable_0.3.6,       data.table_1.18.2.1, hms_1.1.4,,  BiocSingular_1.24.0, ScaledMatrix_1.16.0, xml2_1.5.2,, XVector_0.48.0,     BiocGenerics_0.54.1, pillar_1.11.1,      stringr_1.6.0,      babelgene_22.9,     GSVA_2.2.0,, dplyr_1.2.0,, BiocFileCache_2.16.2, lattice_0.22-7,     bit_4.6.0,,  annotate_1.86.1,    tidyselect_1.2.1,   SingleCellExperiment_1.30.1 Biostrings_2.76.0,  pbapply_1.7-4,      IRanges_2.42.0, SummarizedExperiment_1.38.1 stats4_4.5.2,       Biobase_2.68.0,     statmod_1.5.1,      matrixStats_1.5.0,  pheatmap_1.0.13,    stringi_1.8.7,      UCSC.utils_1.4.0,   codetools_0.2-20,   tibble_3.3.1,       graph_1.86.0,       cli_3.6.5,,  xtable_1.8-8,       Rdpack_2.6.6,       Rcpp_1.1.1,, GenomeInfoDb_1.44.3, dbplyr_2.5.2,      png_0.1-8
library(biomaRt)
library(GSEMA)
library(ActivePathways)
library(limma)
library(writexl)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)

#saving default graphical parameters
opar=par(no.readonly = T)

#Current directory (in rstudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Paper where LUAD data was obtained: https://doi.org/10.1038/nature13385
#Paper of PanCanAtlas: https://doi.org/10.1038/ng.2764
#LinkedOmicsKB link (where data was downloaded from): https://kb.linkedomics.org/download#LUAD

#Read and process omics data ----
exprTumor <- read.delim("LinkedOmicsKB/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt", row.names = 1)
protTumor <- read.delim("LinkedOmicsKB/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt", row.names = 1)
exprNormal <- read.delim("LinkedOmicsKB/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt", row.names = 1)
protNormal <- read.delim("LinkedOmicsKB/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt", row.names = 1)

#Discard genes in pseudoautosomal regions
exprTumor <- exprTumor[grep("PAR", rownames(exprTumor), invert = T),]
exprNormal <- exprNormal[grep("PAR", rownames(exprNormal), invert = T),]

#Remove ENSEMBL version numbers
rownames(exprTumor) <- gsub("\\..*", "", rownames(exprTumor))
rownames(protTumor) <- gsub("\\..*", "", rownames(protTumor))
rownames(exprNormal) <- gsub("\\..*", "", rownames(exprNormal))
rownames(protNormal) <- gsub("\\..*", "", rownames(protNormal))

#Translate ENSEMBL IDs to Gene Symbols
#Query BioMart
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensemble_ids=unique(c(
  rownames(exprTumor),
  rownames(exprNormal),
  rownames(protTumor),
  rownames(protNormal)
))
GP_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                 values=ensemble_ids,mart = mart)
GP_list <- GP_list[!duplicated(GP_list[,1]),]#Remove duplicate names
rownames(GP_list) <- GP_list[,1]

#Function that translates to Gene Symbol and aggregates duplicated genes

ens2sym=function(data,mapping){
  common_ids=intersect(rownames(data),mapping[,2])#get ids from current data matrix
  data=data[rownames(data) %in% mapping[,1],]#remove genes with no mapping
  data$GeneSymbol=mapping[rownames(data),2]#add column with symbol for aggregation
  data_symbol=aggregate(data,list(data$GeneSymbol),median)#aggregate rows with same gene_symbol by median
  data_symbol=na.omit(data_symbol)#remove na
  rownames(data_symbol) <- data_symbol[,1]#add symbol as rownames
  data_symbol=data_symbol[,-1]#Remove the median of genes without symbol
  data_symbol=data_symbol[,-ncol(data_symbol)]#Remove column $GeneSymbol
  return(data_symbol)
}


exprTumor_symbol=ens2sym(exprTumor,GP_list)
exprNormal_symbol=ens2sym(exprNormal,GP_list)
protTumor_symbol=ens2sym(protTumor,GP_list)
protNormal_symbol=ens2sym(protNormal,GP_list)

#Save data for later scripts

save("exprTumor_symbol","exprNormal_symbol","protTumor_symbol","protNormal_symbol",file="Ready/LUAD_expr_prot_21_07.RData")

#Different number of proteins came out (7943 and 7903 for tumor and normal respectively), unify protein matrices with the following vector of common proteins
common_tumor_prot=intersect(rownames(protTumor_symbol),rownames(protNormal_symbol))

#Add group naming to samples
colnames(exprNormal_symbol)=paste("Ctrl",colnames(exprNormal_symbol),sep="_")
colnames(exprTumor_symbol)=paste("Tmr",colnames(exprTumor_symbol),sep="_")
colnames(protNormal_symbol)=paste("Ctrl",colnames(protNormal_symbol),sep="_")
colnames(protTumor_symbol)=paste("Tumor",colnames(protTumor_symbol),sep="_")

#Join groups inside omics
expr_symbol=cbind(exprNormal_symbol,exprTumor_symbol)
prot_symbol=cbind(protNormal_symbol[common_tumor_prot,],protTumor_symbol[common_tumor_prot,])

#Remove genes with all expression values 0 (this only occurs in genes)
expr_symbol=expr_symbol[-which(rowSums(expr_symbol)==0),]

#Clean memory
rm(exprNormal_symbol,exprTumor_symbol,protNormal_symbol,protTumor_symbol,
   exprNormal,exprTumor,protNormal,protTumor,
   GP_list,mart,ensemble_ids,common_tumor_prot,ens2sym)

#Load c2 subcollections (CGP and CP)
msigdb_C2_CGP=msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
msigdb_C2_CP=msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP")
msigdb_C2_sub=rbind(msigdb_C2_CGP,msigdb_C2_CP)
save(msigdb_C2_sub,file="Ready/msigdbr_C2_sub.RData")
#Format genesets as a list
msigdb_C2_sub=split(x=msigdb_C2_sub$gene_symbol,f=msigdb_C2_sub$gs_name)
#remove duplicated gene symbols
msigdb_C2_sub=lapply(msigdb_C2_sub,unique)


#GSEMA ----
#Prepare the labels to distinguish groups
expr_pheno=prot_pheno=data.frame(Condition=c(rep("Healthy",101),rep("Case",110)))
rownames(expr_pheno)=colnames(expr_symbol)
rownames(prot_pheno)=colnames(prot_symbol)
#Vectors to identify groups and labels in each omic
phenoGroups <- c("Condition","Condition")
phenoControls <- c("Healthy", "Healthy")
phenoCases <- c("Case", "Case")
#Perform GSEMA
GSEMA_scores=createObjectMApath(listEX=list(expr_symbol,prot_symbol),
                               listPheno=list(expr_pheno,prot_pheno),
                               namePheno=phenoGroups,
                               expGroups=phenoCases,
                               refGroups=phenoControls,
                               geneSets=msigdb_C2_sub,
                               pathMethod="GSVA",
                               minSize=2,
                               n.cores = 10,#ELIMINAR
                               internal.n.cores = 10)#ELIMINAR

GSEMA_results <- metaAnalysisESpath(objectMApath = GSEMA_scores,
                                   measure = "limma", typeMethod = "REM",
                                   numData = length(GSEMA_scores))

# DPM ----
#Prepare desing matrix for limma
group <- factor(expr_pheno$Condition,levels=c("Healthy","Case"))#Same for protein and expression
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "HealthyVsCase")
#Do limma differential expression
expr_fit=lmFit(expr_symbol,design);expr_fit=eBayes(expr_fit)
prot_fit=lmFit(prot_symbol,design);prot_fit=eBayes(prot_fit)

expr_limma_gene_res=topTable(expr_fit,coef="HealthyVsCase",number=Inf)
prot_limma_gene_res=topTable(prot_fit,coef="HealthyVsCase",number=Inf)
#Obtain results that are common to both omics
dpm_expr_prot_gen_intersect=intersect(rownames(expr_limma_gene_res),rownames(prot_limma_gene_res))
#Store pvalues (nominal, as in their paper)
dpm_pval_matrix=data.frame(row.names = dpm_expr_prot_gen_intersect,
                           rna=expr_limma_gene_res[dpm_expr_prot_gen_intersect,]$P.Value,
                           protein=prot_limma_gene_res[dpm_expr_prot_gen_intersect,]$P.Value)
dpm_pval_matrix=as.matrix(dpm_pval_matrix)
#Store directionality (based on t-values)
dpm_dir_matrix=data.frame(row.names = dpm_expr_prot_gen_intersect,
                          rna=expr_limma_gene_res[dpm_expr_prot_gen_intersect,]$t,
                          protein=prot_limma_gene_res[dpm_expr_prot_gen_intersect,]$t)
dpm_t_matrix=as.matrix(dpm_dir_matrix)#save for later comparisons
dpm_dir_matrix=sign(dpm_dir_matrix)
dpm_dir_matrix=as.matrix(dpm_dir_matrix)
#Store the expected directionality
dpm_constraints_vector=c(1,1)

#This code is to merge p values at gene level
directional_merged_pvals <- merge_p_values(dpm_pval_matrix,
                                           method = "DPM",
                                           dpm_dir_matrix,
                                           dpm_constraints_vector)

#DPM requires pathways to be given in GMT format, which is basically the same as the format given above, but each element of the list has the name of the pathway (in id, no as the name of the element)
#a description and the genes has to be written to a file in a format similar to csv

gmt_msigdb_C2_sub=lapply(names(msigdb_C2_sub), function(x){
  c(x,"NA",msigdb_C2_sub[[x]])
})
writeLines(
  sapply(gmt_msigdb_C2_sub, paste, collapse = "\t"),
  con = "Ready/msigdbr_C2_sub.gmt"
)


dpm_results <- ActivePathways(
  dpm_pval_matrix, gmt = "Ready/msigdbr_C2_sub.gmt",significant = 1,correction_method = "fdr", merge_method = "DPM",
  scores_direction = dpm_dir_matrix, constraints_vector = dpm_constraints_vector, geneset_filter = c(2,3000))#The biggest set in this subcollection has size 2413

save(GSEMA_results,dpm_results,expr_limma_gene_res,prot_limma_gene_res,file="Ready/LUAD_methods_results.RData")

#Results export (tables and graphs)----

#Heatmaps
#Function that takes omics matrices, t-values, p-values, groups, a pathway, color points (in quantile) and file name to produce its heatmap

heatmap_expr_prot=function(expr,prot,
                           expr_t,
                           expr_p,prot_p,groups,pathway,col_quantile,filename){
  #Prepare t-values and extrat the ordering
  expr_t=expr_t[match(pathway,rownames(expr_t)),,F]
  expr_t_p_order=order(expr_t[,1],decreasing = T)
  expr_t=expr_t[expr_t_p_order,,F]
  #Prepare and order p-values
  expr_p=expr_p[match(pathway,rownames(expr_p)),,F]
  expr_p=expr_p[expr_t_p_order,,F]
  #Prepare and order protein p-values
  prot_p=prot_p[match(pathway,rownames(prot_p)),,F]
  rownames(prot_p)=pathway
  prot_p=prot_p[match(rownames(expr_p),rownames(prot_p)),,F]
  #Prepare significance coding
  expr_sig=symnum(expr_p$adj.P.Val,c(0,0.001,0.01,0.05,0.1,1),c("***","**","*","·"," "),cor=F,na=" ")
  prot_sig=symnum(prot_p$adj.P.Val,c(0,0.001,0.01,0.05,0.1,1),c("***","**","*","·"," "),cor=F,na=" ")
  #Prepare omics matrices with only the expression of the path and the ordering of expression t-values
  expr=expr[match(pathway,rownames(expr)),]
  expr=as.matrix(expr[na.omit(match(rownames(expr_t),rownames(expr))),])

  prot=prot[match(pathway,rownames(prot)),]
  prot=as.matrix(prot[match(rownames(expr_t),rownames(prot)),])
  rownames(prot)=rownames(expr_t)
  #Object to plot t-value vertical barplot
  t_val_annotation <- rowAnnotation(
    `t-value` = anno_barplot(expr_t, gp = gpar(fill = "darkgrey")),
    annotation_name_side = "top"
  )
  #Object to plot sampling groups
  group_annotation <- HeatmapAnnotation(
    Cohort = groups,
    col = list(Cohort = c("Healthy" = "royalblue", "Case" = "tan1")),
    show_annotation_name = FALSE
  )
  #Coloring
  expr_vect=as.numeric(t(scale(t(expr))))
  prot_vect=as.numeric(t(scale(t(expr))))
  expr_col_points=c(quantile(expr_vect,col_quantile),0,quantile(expr_vect,1-col_quantile))
  prot_col_points=c(quantile(prot_vect,col_quantile,na.rm=T),0,quantile(prot_vect,1-col_quantile,na.rm=T))
  col_fun1 <- colorRamp2(expr_col_points, c("blue", "white", "red"))
  col_fun2 <- colorRamp2(prot_col_points, c("blue", "white", "red"))
  #object to produce Expression heatmap
  ht_expr <- Heatmap(t(scale(t(expr))),
                     name = "Expression",
                     col = col_fun1,
                     heatmap_legend_param = list(
                       at = expr_col_points,
                       labels = as.character(round(expr_col_points))),
                     cluster_rows = FALSE,
                     row_labels = as.character(expr_sig),
                     cluster_columns = F,
                     column_title = "Gene Expression",
                     show_column_names = F,
                     show_row_names = T,
                     top_annotation = group_annotation)

  expr_stars_anno <- rowAnnotation(
    Expr_Sig = anno_text(as.character(expr_sig),
                         gp = gpar(fontface = "bold"),
                         just = "center",
                         location = 0.5)
  )

  ht_prot <- Heatmap(t(scale(t(prot))),
                     name = "Proteomics",
                     col = col_fun2,
                     heatmap_legend_param = list(
                       at = prot_col_points,
                       labels = as.character(round(prot_col_points))),
                     cluster_rows = FALSE,
                     row_labels = as.character(prot_sig),
                     cluster_columns = F,
                     column_title = "Proteomics",
                     show_column_names = F,
                     na_col = "grey90",
                     top_annotation = group_annotation)
  pdf(filename, width = 13, height = 9,bg="white")
  #png(file=filename, units="in", width=11, height=8.5, res=1200)
  draw(ht_expr+expr_stars_anno+t_val_annotation+ht_prot, row_title = "Pathway Genes (Ordered by t-value)")
  dev.off()
}
#Produce heatmaps
pathways_4_heatmaps=paste0(c("1","s11","s12","s13","s14","s15","2","3","s31","s32","s33"),"_",
                          c("YU_BAP1_TARGETS",
                            "MEISSNER_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3",
                            "MIKKELSEN_ES_HCP_WITH_H3K27ME3",
                            "TURJANSKI_MAPK14_TARGETS",
                            "MEBARKI_HCC_PROGENITOR_WNT_DN_CTNNB1_DEPENDENT_BLOCKED_BY_FZD8CRD",
                            "BENPORATH_ES_2",
                            "GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_BLACK_DN",
                            "GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_CYAN_UP",
                            "VISALA_AGING_LYMPHOCYTE_DN",
                            "NABA_MATRISOME_POORLY_METASTATIC_MELANOMA_TUMOR_CELL_DERIVED",
                            "WATANABE_ULCERATIVE_COLITIS_WITH_CANCER_UP"))

lapply(pathways_4_heatmaps, function(x){
  heatmap_expr_prot(expr_symbol,prot_symbol,
                    expr_limma_gene_res[,"t",F],
                    expr_limma_gene_res[,"adj.P.Val",F],prot_limma_gene_res[,"adj.P.Val",F],
                    expr_pheno$Condition,
                    getElement(msigdb_C2_sub,sub("^[^_]*_", "", x)),
                    0.001,paste0("Figures/Heatmaps/",x,".pdf"))
})



#Save genes and proteins present
genes_in_expr=rownames(expr_limma_gene_res)
genes_in_prot=rownames(prot_limma_gene_res)
#Save their t values
expr_limma_gene_res_t=expr_limma_gene_res$t
names(expr_limma_gene_res_t)=rownames(expr_limma_gene_res)
prot_limma_gene_res_t=prot_limma_gene_res$t
names(prot_limma_gene_res_t)=rownames(prot_limma_gene_res)
#take what is needed from dpm (term id and adjusted p-value)
df_dpm_results=as.data.frame(dpm_results[,c(1,3)])
rownames(df_dpm_results)=dpm_results$term_id

#Now we are groing to get pathways that lose information in p-value meta-analysis because of genes removed by not having matching protein.
#Get which genes do have a corresponding protein
ind_genes_with_prot=rownames(expr_limma_gene_res) %in% rownames(prot_limma_gene_res)
names(ind_genes_with_prot)=rownames(expr_limma_gene_res)

#Get, for each pathway, the contribution of all the genes that are missing (calculated as the sum of t values of all genes with no protein with respect to all genes)
path_without_sig_genes_by_prot=t(sapply(msigdb_C2_sub, function(x){
  ind_genes_in_set=ind_genes_with_prot[x]
  ind_genes_in_set=ind_genes_in_set[!is.na(ind_genes_in_set)]
  ind_genes_in_set_contribution=abs(expr_limma_gene_res[names(ind_genes_in_set),"t"])
  names(ind_genes_in_set_contribution)=names(names(ind_genes_in_set))
  rel_contribution=sum(ind_genes_in_set_contribution[which(ind_genes_in_set==0)])/sum(ind_genes_in_set_contribution)
  c(rel_contribution,length(ind_genes_in_set))
}))

#Join the contribution (and number of genes in the set) with the p-values of GSEMA and DPM
path_without_sig_genes_by_prot_joined=data.frame(rownames(path_without_sig_genes_by_prot),
                                                 path_without_sig_genes_by_prot,
                                                 df_dpm_results[rownames(path_without_sig_genes_by_prot),"adjusted_p_val"],
                                                 GSEMA_results[rownames(path_without_sig_genes_by_prot),"FDR"]
)


colnames(path_without_sig_genes_by_prot_joined)=c("Path","Missing_Contribution","Path_size","DPM_FDR","GSEMA_FDR")
#Plot results

color_size=colorRampPalette(c('lightblue', 'blue'))
temp_ind=as.numeric(cut(log(path_without_sig_genes_by_prot_joined[,3]),breaks = 100))#index for color ordered by size
color_size=color_size(100)[temp_ind]

pdf("Figures/Jurado-Bascon_Fig5.pdf", width = 13, height = 9,bg="white")
#png(file="Figures/Jurado-Bascon_Fig5.png", units="in", width=11, height=8.5, res=1200)
par(mfrow=c(2,1),mar=c(4,4,0,0)+0.5,xpd=T)
plot(1-path_without_sig_genes_by_prot_joined[,2],
     path_without_sig_genes_by_prot_joined[,4],
     col=color_size,
     xlab="Cross-omics feature overlap contribution",ylab="DPM FDR",pch=19,cex=0.5)
legend("topright",inset=c(0.05,0),legend=c("Size~10","Size~1700"),pch=19,col=color_size[c(which.min(temp_ind),which.max(temp_ind))],bty="n")
plot(1-path_without_sig_genes_by_prot_joined[,2],
     path_without_sig_genes_by_prot_joined[,5],
     col=color_size,
     xlab="Cross-omics feature overlap contribution",ylab="GSEMA FDR",pch=19,cex=0.5)
par(opar)
dev.off()

#Get pathways that have genes and their proteins following the same directionality (concordance)
#As DPM only takes into account genes that have associated proteins, only those genes will be considered
matched_expr_limma_gene_res=expr_limma_gene_res[rownames(prot_limma_gene_res),]
#Get a named vector of the t-values for genes and proteins
matched_expr_t=setNames(matched_expr_limma_gene_res$t,rownames(matched_expr_limma_gene_res))
prot_t=setNames(prot_limma_gene_res$t,rownames(prot_limma_gene_res))
#Get a named vector of the sign for genes and proteins
matched_expr_sign=matched_expr_t>0
prot_sign=prot_t>0
#Get the concordance, measured as the sum of absolute t-values of concordance pair divided by the sum of all t-value
concordance_proportion=(sapply(msigdb_C2_sub, function(x){
  x=x[x %in% names(matched_expr_sign)]#Remove missing values
  matched_expr_sign_current=matched_expr_sign[x]#get sign for expr
  prot_sign_current=prot_sign[x]#get sign for prot
  abs_t_current=abs(cbind(matched_expr_t[x],prot_t[x]))#get absolute expression from both
  ind_concordance_terms=(matched_expr_sign_current+prot_sign_current)!=1#get which genes and proteins agree on sigh
  sum((rowSums(abs_t_current))[ind_concordance_terms])/sum(rowSums(abs_t_current))#Return the contribution of genes with equal sign
}))

concordance_proportion_joined=data.frame(names(concordance_proportion),
                                         concordance_proportion,
                                         df_dpm_results[names(concordance_proportion),"adjusted_p_val"],
                                         GSEMA_results[names(concordance_proportion),"FDR"]
)
colnames(concordance_proportion_joined)=c("Path","Concordance","DPM_FDR","GSEMA_FDR")

pdf("Figures/Jurado-Bascon_Fig7.pdf", width = 13, height = 9,bg="white")
#png(file="Figures/Jurado-Bascon_F7.png", units="in", width=11, height=8.5, res=1200)
par(mfrow=c(2,1),mar=c(4,4,0,0)+0.5)
plot(concordance_proportion_joined[,2],
     concordance_proportion_joined[,3],
     col=color_size,
     xlab="Cross-omics directional concordance",ylab="DPM FDR",pch=19,cex=0.5)
legend("topleft",inset=c(0.05,0),legend=c("Size~10","Size~1700"),pch=19,col=color_size[c(which.min(temp_ind),which.max(temp_ind))],bty="n")
plot(concordance_proportion_joined[,2],
     concordance_proportion_joined[,4],
     col=color_size,
     xlab="Cross-omics directional concordance",ylab="GSEMA FDR",pch=19,cex=0.5)
par(opar)
dev.off()
#Get pathways that disagree in directionality (they may appear as relevant by their individual t values and significance, but the sign of the t-values does not show a common directionality for the set)
enrichment_proportion=(sapply(msigdb_C2_sub, function(x){
  x=x[x %in% names(matched_expr_sign)]#Remove missing values
  t_current=cbind(matched_expr_t[x],prot_t[x])#get expression from both
  sum(rowSums(t_current))/sum(rowSums(abs(t_current)))#Return the directionality of the pathway
}))

enrichment_proportion_joined=data.frame(names(enrichment_proportion),
                                        enrichment_proportion,
                                        df_dpm_results[names(enrichment_proportion),"adjusted_p_val"],
                                        GSEMA_results[names(enrichment_proportion),"FDR"]
)

colnames(enrichment_proportion_joined)=c("Path","Concordance","DPM_FDR","GSEMA_FDR")

pdf("Figures/Jurado-Bascon_Fig8.pdf", width = 13, height = 9,bg="white")
#png(file="Figures/F3_Disagreement_Contribution.png", units="in", width=11, height=8.5, res=1200)
par(mfrow=c(2,1),mar=c(4,4,0,0)+0.5)
plot(enrichment_proportion_joined[,2],
     enrichment_proportion_joined[,3],
     col=color_size,
     xlab="Pathway regulatory homogeneity",ylab="DPM FDR",pch=19,cex=0.5)
legend("topleft",inset=c(0.05,0),legend=c("Size~10","Size~1700"),pch=19,col=color_size[c(which.min(temp_ind),which.max(temp_ind))],bty="n")
plot(enrichment_proportion_joined[,2],
     enrichment_proportion_joined[,4],
     col=color_size,
     xlab="Pathway regulatory homogeneity",ylab="GSEMA FDR",pch=19,cex=0.5)
par(opar)
dev.off()

#Results table
#Obtain the posible directionality of the sets, in this case, the directionality of its significant elements, mean of t-value and so on (for results exploration)
sig_genes_in_expr=rownames(expr_limma_gene_res[expr_limma_gene_res$P.Value<0.05,])
sig_genes_in_prot=rownames(prot_limma_gene_res[prot_limma_gene_res$P.Value<0.05,])

expr_limma_sig_gene_t=expr_limma_gene_res$t[expr_limma_gene_res$P.Value<0.05]
names(expr_limma_sig_gene_t)=sig_genes_in_expr
prot_limma_sig_gene_t=prot_limma_gene_res$t[prot_limma_gene_res$P.Value<0.05]
names(prot_limma_sig_gene_t)=rownames(prot_limma_gene_res[prot_limma_gene_res$P.Value<0.05,])

expr_path_dir=t(sapply(msigdb_C2_sub, function(x){#Proportion of significant t-values that are over expressed in genes
  temp=expr_limma_sig_gene_t[x]
  temp=temp[!is.na(temp)]
  c(sum(temp>0)/(length(temp)),mean(temp))
}))

prot_path_dir=t(sapply(msigdb_C2_sub, function(x){#Proportion of significant t-values that are over expressed in proteins
  temp=prot_limma_sig_gene_t[x]
  temp=temp[!is.na(temp)]
  c(sum(temp>0)/(length(temp)),mean(temp))
}))

#Get other relevant information from pathways
genes_present_in_gene_sets=t(sapply(msigdb_C2_sub, function(x){#Proportion of genes present in the pathway
  c(sum(x %in% genes_in_expr)/length(x),sum(x %in% sig_genes_in_expr)/length(x))
}))

prots_present_in_gene_sets=t(sapply(msigdb_C2_sub, function(x){#Proportion of proteins present in the pathway
  c(sum(x %in% genes_in_prot)/length(x),sum(x %in% sig_genes_in_prot)/length(x))
}))

path_size=sapply(msigdb_C2_sub, length)#Size of the pathway

genes_up_in_gene_sets=sapply(msigdb_C2_sub, function(x){#Most relevant (as greater t value) genes in the pathway
  temp=sort(expr_limma_gene_res_t[x],T)[1:10]
  paste(names(temp), temp, collapse = " ; ")
})
#The same for lower t-value an for proteins
genes_down_in_gene_sets=sapply(msigdb_C2_sub, function(x){
  temp=sort(expr_limma_gene_res_t[x],F)[1:10]
  paste(names(temp), temp, collapse = " ; ")
})

prots_up_in_gene_sets=sapply(msigdb_C2_sub, function(x){
  temp=sort(prot_limma_gene_res_t[x],T)[1:10]
  paste(names(temp), temp, collapse = " ; ")
})

prots_down_in_gene_sets=sapply(msigdb_C2_sub, function(x){
  temp=sort(prot_limma_gene_res_t[x],F)[1:10]
  paste(names(temp), temp, collapse = " ; ")
})

enrichment_proportion_joined$Genes.up=genes_up_in_gene_sets[enrichment_proportion_joined$Path]
enrichment_proportion_joined$Genes.down=genes_down_in_gene_sets[enrichment_proportion_joined$Path]
enrichment_proportion_joined$Prots.up=prots_up_in_gene_sets[enrichment_proportion_joined$Path]
enrichment_proportion_joined$Prots.down=prots_down_in_gene_sets[enrichment_proportion_joined$Path]

pathway_summary=data.frame(path_size,
                           genes_up_in_gene_sets,genes_down_in_gene_sets,genes_present_in_gene_sets,expr_path_dir,
                           prots_up_in_gene_sets,prots_down_in_gene_sets,prots_present_in_gene_sets,prot_path_dir)
colnames(pathway_summary)=c("Size",
                            "Genes.up","Genes.dn","%genes", "%sig.genes", "%OE.genes", "mean.t.genes",
                            "Prots.up","Prots.dn","%prots", "%sig.prots", "%OE.prots", "mean.t.prots")
paths=rownames(pathway_summary)
#Join results in excel
all_methods_results=data.frame(
  Path=paths,
  pathway_summary,
  Missing_Contribution=path_without_sig_genes_by_prot_joined[,2],
  Concordance_Contribution=concordance_proportion_joined[,2],
  Directionality_agreement=enrichment_proportion_joined[,2],
  DPM_Path=df_dpm_results[paths,"term_id"],
  DPM_FDR=df_dpm_results[paths,"adjusted_p_val"],
  GSEMA_Path=GSEMA_results[paths,"Pathway"],
  GSEMA_Eff=GSEMA_results[paths,"Com.ES"],
  GSEMA_FDR=GSEMA_results[paths,"FDR"]
)

all_methods_results_sheets=list(
  all_methods_results=all_methods_results,
  missing_contribution=path_without_sig_genes_by_prot_joined,
  disagreement=enrichment_proportion_joined
)

#write_xlsx(all_methods_results_sheets,"Tables/LUAD_Method_Compare_C2_sub_22_07_26.xlsx")
