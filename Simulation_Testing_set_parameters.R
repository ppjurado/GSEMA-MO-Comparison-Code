#Loading packages
library(msigdbr)
library(doParallel)
library(doRNG)
library(GSEMA)
library(ActivePathways)

#saving default graphical parameters
opar=par(no.readonly = T)

#Current directory (in rstudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Save functions and parameters used for simulation: ----
n_iter=100
n_cores=10

#Function that executes GSEMA inside simulation
iteration_GSEMA=function(data,group,sets){
  #Prepare the labels to distinguish groups
  prot_pheno=expr_pheno=data.frame(Condition=group)
  rownames(expr_pheno)=colnames(data[[1]])
  rownames(prot_pheno)=colnames(data[[2]])
  #Vectors to identify groups and labels
  phenoGroups <- c("Condition","Condition")
  phenoControls <- c("Healthy", "Healthy")
  phenoCases <- c("Case", "Case")
  #obtain scores for GSEMA
  GSEMA_scores=createObjectMApath(listEX=data,
                                 listPheno=list(expr_pheno,prot_pheno),
                                 namePheno=phenoGroups,
                                 expGroups=phenoCases,
                                 refGroups=phenoControls,
                                 geneSets=sets,
                                 pathMethod="GSVA")
  GSEMA_results <- metaAnalysisESpath(objectMApath = GSEMA_scores,
                                     measure = "limma", typeMethod = "REM", missAllow = 0.3,
                                     numData = length(GSEMA_scores))
  return(GSEMA_results[,c(2,3,7)])#Return effect size, its variance and nominal p-value
}
#Function that executes DPM inside simulation
iteration_DPM=function(data,group,sets){
  #Prepare desing matrix for limma
  group=factor(group,levels=c("Healthy","Case"))
  design <- model.matrix(~ group)
  colnames(design) <- c("Intercept", "HealthyVsCase")
  #Do limma differential expression
  expr_fit=lmFit(data[[1]],design);expr_fit=eBayes(expr_fit)
  prot_fit=lmFit(data[[2]],design);prot_fit=eBayes(prot_fit)
  expr_limma_gene_res=topTable(expr_fit,coef="HealthyVsCase",number=Inf)
  prot_limma_gene_res=topTable(prot_fit,coef="HealthyVsCase",number=Inf)

  #Obtain results that are common to both omics (as needed by DPM)
  dpm_expr_prot_gen_intersect=intersect(rownames(expr_limma_gene_res),rownames(prot_limma_gene_res))
  #Store pvalues
  dpm_pval_matrix=data.frame(row.names = dpm_expr_prot_gen_intersect,
                             rna=expr_limma_gene_res[match(dpm_expr_prot_gen_intersect,rownames(expr_limma_gene_res)),]$P.Value,
                             protein=prot_limma_gene_res[match(dpm_expr_prot_gen_intersect,rownames(prot_limma_gene_res)),]$P.Value)
  dpm_pval_matrix=as.matrix(dpm_pval_matrix)
  #Store directionality (based on t-values)
  dpm_dir_matrix=data.frame(row.names = dpm_expr_prot_gen_intersect,
                            rna=expr_limma_gene_res[match(dpm_expr_prot_gen_intersect,rownames(expr_limma_gene_res)),]$t,
                            protein=prot_limma_gene_res[match(dpm_expr_prot_gen_intersect,rownames(prot_limma_gene_res)),]$t)
  dpm_dir_matrix=sign(dpm_dir_matrix)
  dpm_dir_matrix=as.matrix(dpm_dir_matrix)
  #Store the expected directionality agreement between genes and proteins
  dpm_constraints_vector=c(1,1)
  dpm_pval_matrix[is.na(dpm_pval_matrix)]=1#Some iterations produce NA in limma, this happens when doing a lot of simulations with effects, not in the normal analysis or in false positive simulation
  dpm_dir_matrix[is.na(dpm_dir_matrix)]=0
  dpm_results <- ActivePathways(
    dpm_pval_matrix, gmt = "Ready/simulation_sets.gmt",significant = 1,correction_method = "none", merge_method = "DPM",
    scores_direction = dpm_dir_matrix, constraints_vector = dpm_constraints_vector, geneset_filter = c(5, 10000))
  return(as.data.frame(dpm_results[,c(1,3)]))#Return path name and nominal p-value (parameter correction_method = "none")
}
#Function that executes one iteration of a simulation (give data generated with concrete parameters, shuffled group labels and sets to study enrichment)
iteration=function(data,group,sets){#Join the results
  dpm_res=iteration_DPM(data,group,sets)
  GSEMA_res=iteration_GSEMA(data,group,sets)
  it_res=cbind(dpm_res,GSEMA_res)
  colnames(it_res)=c("Path","DPM_pv","GSEMA_eff","GSEMA_eff_var","GSEMA_pv")
  return(it_res)
}
#Function that, given the expression matrices, applies the specified effect to randomly selected proportion of genes and their proteins, with a determined proportion of positive and negative signs
apply_effect=function(data,group,sets,effect,prop_s,prop_b,prop_d=1){
  #Prepare easy access for data
  expr=data[[1]]
  prot=data[[2]]
  #Sample the genes and proteins in the sets affected by the positive effect
  c_g_sets_pos=lapply(sets, function(x){
    ns_pos=round(length(x)*prop_s*prop_d)
    sample(x,ns_pos)})
  c_g_sets_neg=lapply(names(sets), function(x){
    ns_neg=round(length(sets[[x]])*prop_s*(1-prop_d))
    y=sets[[x]][!(sets[[x]]%in%c_g_sets_pos[[x]])]
    sample(y,ns_neg)})
  c_p_sets_pos=lapply(c_g_sets_pos, function(x){rownames(prot)[rownames(prot) %in% x]})#Get the proteins associated with the sampled genes with positive effect
  c_p_sets_neg=lapply(c_g_sets_neg, function(x){rownames(prot)[rownames(prot) %in% x]})#Get the proteins associated with the sampled genes with negative effect
  #Number of genes in background affected by positive and negative effects
  nb_pos=round(length(genes_not_in_sets)*prop_b*prop_d)
  nb_neg=round(length(genes_not_in_sets)*prop_b*(1-prop_d))
  #sample genes and proteins in background
  c_g_back_pos=sample(genes_not_in_sets,nb_pos)#Sample the genes in the background affected by the positive effect
  c_g_back_neg=sample(genes_not_in_sets[!(genes_not_in_sets%in%c_g_back_pos)],nb_neg)#Extract the genes with positive effect, and from those sample the assigned to negative effect
  c_p_back_pos=prots_not_in_sets[prots_not_in_sets %in% c_g_back_pos]#get the proteins associated with the sampled background genes
  c_p_back_neg=prots_not_in_sets[prots_not_in_sets %in% c_g_back_neg]
  cases=group=="Case"

  #First, positive effect
  ind=match(unique(c(unlist(c_g_sets_pos),c_g_back_pos)),rownames(expr))
  expr[ind,cases]=expr[ind,cases]+#cases in those genes will be them summed by
    effect*apply(expr[ind,], 1, sd,na.rm=T)#the effect times their standard deviation
  #The same with protein
  ind=match(unique(c(unlist(c_p_sets_pos),c_p_back_pos)),rownames(prot))
  prot[ind,cases]=prot[ind,cases]+
    effect*apply(prot[ind,], 1, sd,na.rm=T)
  #Now negative
  ind=match(unique(c(unlist(c_g_sets_neg),c_g_back_neg)),rownames(expr))
  expr[ind,cases]=expr[ind,cases]-#cases in those genes will be them summed by
    effect*apply(expr[ind,], 1, sd,na.rm=T)#the effect times their standard deviation
  #The same with protein
  ind=match(unique(c(unlist(c_p_sets_neg),c_p_back_neg)),rownames(prot))
  prot[ind,cases]=prot[ind,cases]-
    effect*apply(prot[ind,], 1, sd,na.rm=T)
  return(list(expr,prot))
}
#Function that parallelises the simulation. It takes a combination of effect parameters, applies the effect and shuffles the labels and saves relevant results from DPM and GSEMA
scenario=function(data,group,sets,n_iter,n_cores,effect,prop_s,prop_b,prop_d=1){
  t1=Sys.time()#Count time of start
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)#Prepare clusters, for windows
  on.exit(stopCluster(cl))#ensure execution of stopCluster
  sim_res=foreach(i=1:n_iter,#vector to iterate over
                  .inorder=F,#tells foreach if the order of the iterations matter (F, as it takes less time and it does not matter)
                  .packages=c("GSEMA","limma","ActivePathways"),#Needed packages
                  .export=c("iteration","iteration_DPM","iteration_GSEMA","apply_effect","genes_not_in_sets","prots_not_in_sets"),#workers by default do not know about earlier defined functions and objects, this ensures they know
                  .final=function(x){#Function to join the results of the iterations
                    it_res=do.call(rbind,x)#x is the result of the chunk after in %dorng%
                    return(split(it_res[,-1],it_res$Path))#Separate results by set
                  }) %dorng% {#Even when doing parallel computing, operator %dorng% ensures each cluster gets a different chain of random numbers
                    group=sample(group)#randomize groups
                    data_eff=apply_effect(data,group,sets,effect,prop_s,prop_b,prop_d)#Apply effect
                    return(iteration(data_eff,group,sets))#execute iteration
                  }
  attr(sim_res,"rng")=NULL#Remove seeds provided by doRNG
  attr(sim_res,"doRNG_version")=NULL#Remove doRNG version
  t2=Sys.time()#Count time of finish
  message("Execution time: ",round(t2-t1,2),"\n number of interations:",n_iter,"\n Parameters: eff=",round(effect,3),", prop_s=",prop_s,", prop_b=",prop_b,", prop_d=",prop_d,"\n")#Time of execution and parameters of scenario
  return(sim_res)
}

#Load data and sets ----

#Load msigdb C2 whole genesets database (as we will perform a search based on the characteristics of pathway, without interpreting them, we only look for which set profiles are prone to false positives or lack of statistical power)

msigdb_C2=msigdbr(species = "Homo sapiens", collection = "C2")
save(msigdb_C2,file="Ready/msigdbr_C2.RData")
msigdb_C2=split(x=msigdb_C2$gene_symbol,f=msigdb_C2$gs_name)
#remove duplicated gene symbols
msigdb_C2=lapply(msigdb_C2,unique)

#Load LUAD data (saved in LUAD_Enrichment_MO_Analisys)
load("Ready/LUAD_expr_prot.RData")

#Simulations will be performed only with controls, as they have enough samples + the matrix for simulation will be more homogeneus
expr_symbol=exprNormal_symbol[-which(rowSums(exprNormal_symbol)==0),]
prot_symbol=protNormal_symbol
rm(exprTumor_symbol,protTumor_symbol,exprNormal_symbol,protNormal_symbol)

#Calculating set info ----

#We remove gene sets that have less than 95% of genes present (that is only 179 out of 7561)
proportion_genes_present=sapply(msigdb_C2, function(x){
  sum(x %in% rownames(expr_symbol))/length(x)
})
msigdb_C2=msigdb_C2[-which(proportion_genes_present<0.95)]
#For ease of work in the simulation, we remove from all sets the genes that do not appear in the expression matrix
msigdb_C2=lapply(msigdb_C2, function(x){
  x[x %in% rownames(expr_symbol)]
})


#Get the proportion of proteins present
genes_fully_present=intersect(rownames(expr_symbol),rownames(prot_symbol))
proportion_prots_present=sapply(msigdb_C2, function(x){
  sum(x %in% genes_fully_present)/length(x)
})
#Get the mean absolute correlation inside the sets for genes and proteins
corr_prots_genes=t(sapply(msigdb_C2, function(x){
  ind_expr=match(x,rownames(expr_symbol))
  ind_prot=match(x,rownames(prot_symbol))
  c_cor_genes=cor(t(expr_symbol[ind_expr,]))
  c_cor_prots=cor(t(prot_symbol[ind_prot,]))
  c(mean(abs(c_cor_genes[upper.tri(c_cor_genes)]),na.rm=T),mean(abs(c_cor_prots[upper.tri(c_cor_prots)]),na.rm=T))
}))
#Get the mean standard deviation inside the sets for genes and proteins
sds_prots_genes=t(sapply(msigdb_C2, function(x){
  ind_expr=match(x,rownames(expr_symbol))
  ind_prot=match(x,rownames(prot_symbol))
  c(
    mean(apply(expr_symbol[ind_expr,],1,sd,na.rm=T),na.rm=T),
    mean(apply(prot_symbol[ind_prot,],1,sd,na.rm=T),na.rm=T)
  )
}))
#Get the size of the sets
sets_size=sapply(msigdb_C2, length)
#Join information of sets
set_info=data.frame(sets_size,proportion_prots_present,corr_prots_genes,sds_prots_genes)
names(set_info)=c("size","prot.prop","Expr.Cor","Prot.Cor","Expr.Sd","Prot.Sd")
#We are going to filter extreme cases, and go with sets that comply with the following restrictions
#sets with size smaller than 10 (83%)
ind_small=set_info$size>10
#For simplicity, we are going to work with sets that have similar inner correlation in genes and proteins, that is: (58%)
ind_sim_cor=abs(set_info$Expr.Cor-set_info$Prot.Cor)<0.05
ind_sim_cor[is.na(ind_sim_cor)]=F
#the same with sets that have similar standard deviations in genes and proteins, that is: (50%)
ind_sim_sd=abs(set_info$Expr.Sd-set_info$Prot.Sd)<0.2
ind_sim_sd[is.na(ind_sim_sd)]=F
#minimun presence of proteins, a proportion of 0.1 (97%)
ind_small_prot_prop=set_info$prot.prop>0.10
#Combine all filters (26%, that is, 1917 sets)
ind_filt=ind_small & ind_sim_cor & ind_sim_sd & ind_sim_cor

set_info=set_info[ind_filt,]

#TESTING SIZE ----
#First, we are going to test if the size of the set matters, and we will keep the rest of the variables at "common values". Taking the middle 50% of sets in each variable (apart from size) leads to more than 300 sets, therefore we can have a closer interval
#Taking the middle 30% we get 95 sets
cut_set_info=data.frame(
  Size=cut(set_info$size,quantile(set_info$size,seq(0,1,,11)),dig.lab=2,ordered_result = T),
  prot.prop=cut(set_info$prot.prop,quantile(set_info$prot.prop,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  Expr.Cor=cut(set_info$Expr.Cor,quantile(set_info$Expr.Cor,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  Expr.Sd=cut(set_info$Expr.Sd,quantile(set_info$Expr.Sd,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  row.names = rownames(set_info)
)
#common values means that, as the variables are cut by percentiles, the middle 30% of their distribution
#Generate a contingency array to check if we have enough sets for testing
arr_cont=table(cut_set_info)
apply(arr_cont[,2:3,2:3,2:3], 1, sum)
#Now, obtain those sets
mid_levels=t(sapply(cut_set_info[,-1], levels))[,2:3]
mid_cut_set_info=cut_set_info[cut_set_info$prot.prop %in% mid_levels[1,] &
                                cut_set_info$Expr.Cor %in% mid_levels[2,] &
                                cut_set_info$Expr.Sd %in% mid_levels[3,],]
#Return to original size values
mid_cut_set_info$Size=set_info[row.names(mid_cut_set_info),1]
mid_cut_set_info=mid_cut_set_info[order(mid_cut_set_info$Size),]
#Now we have 93 sets, with common values at the rest of the variables, ordered by size, lets analyse them
c_sets=msigdb_C2[rownames(mid_cut_set_info)]
#Set sets for DPM
writeLines(
  sapply(lapply(names(c_sets), function(x){c(x,"NA",c_sets[[x]])}), paste, collapse = "\t"),
  con = "Ready/simulation_sets.gmt"
)

#See which genes and proteins are outside the sets
genes_not_in_sets=rownames(expr_symbol)[which(!(rownames(expr_symbol) %in% unlist(c_sets)))]
prots_not_in_sets=rownames(prot_symbol)[which(!(rownames(prot_symbol) %in% unlist(c_sets)))]
#Define color gradient
grad=colorRampPalette(c("lightblue","darkblue"))

#Now the simulation is run with all set parameters to common values except for size, for which we have a good range of parameters
#We are going to apply 3 different effect sizes and save every p.value plot
#In prior versions of this script, seq(0,0.5,,10) values for effect size where considered, for the article, the 3 values bellow where taken as representatives
#save results in list:
sim_size_x_eff=vector("list")
eff=qnorm(0.5+c(0,0.15,0.4)/2)#We take z values that differ from the mean in a given "percentage", for example, qnorm(0.5+0.25/2)) represents a diferentiation of the mean higher than 25% of the population (this z value is after, in function apply_effect) multiplied by the standard deviation of each gene (de-standarizing the scale)
for(i in 1:length(eff)){
  groups=c(rep("Healthy",50),rep("Case",51))#Set groups
  set.seed(2026)#set seed
  sim_size_x_eff[[i]]=scenario(list(expr_symbol,prot_symbol),groups,c_sets,n_iter,n_cores,eff[i],0.5,0.25,0)#scenario outputs a list, each element is a set in c_sets, each containing the results (pd of dpm and GSEMA) for the 100 iterations
}

for(i in 1:length(eff)){
  c_sim=sim_size_x_eff[[i]]
  #In this case we color proportion of proteins present
  col_grad=grad(100)[as.numeric(cut(set_info[match(names(c_sets),rownames(set_info)),"prot.prop"],100))]
  png(paste0("Figures/STsp/size/",i,"_eff",c(0,0.15,0.4)[i],".png"), units="in", width=11, height=8.5, res=300)
  par(mfrow=c(2,1),mar=c(1.5,4,0,0)+0.5,xpd=t)
  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("DPM_pv")]}),#the sapply takes each of the sets and takes, in this case, the p-value for dpm, therefore getting a matrix that contains the sets (ordered by size) in columns, iterations in rows, the matrix cells containing the dpm p-values
          col=col_grad,names=sapply(c_sets,length),ylim=c(0,1),ylab="DPM p-value")#Change name of set for its size
  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("GSEMA_pv")]}),
          col=col_grad,names=sapply(c_sets,length),ylim=c(0,1),ylab="GSEMA p-value")
  par(opar)
  dev.off()
}


names(sim_size_x_eff)=paste0("eff=",round(eff,2))
save(sim_size_x_eff,file="Ready/STsp_sim_size_x_eff.RData")
#Repeat for proportion of proteins in the set, deviation and correlation
#TESTING PROPORTION OF PROTEINS IN SET ----
cut_set_info=data.frame(
  Size=cut(set_info$size,quantile(set_info$size,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  prot.prop=cut(set_info$prot.prop,quantile(set_info$prot.prop,seq(0,1,,11)),dig.lab=2,ordered_result = T),
  Expr.Cor=cut(set_info$Expr.Cor,quantile(set_info$Expr.Cor,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  Expr.Sd=cut(set_info$Expr.Sd,quantile(set_info$Expr.Sd,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  row.names = rownames(set_info)
)
arr_cont=table(cut_set_info)
apply(arr_cont[2:3,,2:3,2:3], 2, sum)
sum(apply(arr_cont[2:3,,2:3,2:3], 2, sum))#71 sets


mid_levels=t(sapply(cut_set_info[,-2], levels))[,2:3]
mid_cut_set_info=cut_set_info[cut_set_info$Size %in% mid_levels[1,] &
                                cut_set_info$Expr.Cor %in% mid_levels[2,] &
                                cut_set_info$Expr.Sd %in% mid_levels[3,],]

mid_cut_set_info$prot.prop=set_info[row.names(mid_cut_set_info),2]
mid_cut_set_info=mid_cut_set_info[order(mid_cut_set_info$prot.prop),]

c_sets=msigdb_C2[rownames(mid_cut_set_info)]

writeLines(
  sapply(lapply(names(c_sets), function(x){c(x,"NA",c_sets[[x]])}), paste, collapse = "\t"),
  con = "Ready/simulation_sets.gmt"
)

genes_not_in_sets=rownames(expr_symbol)[which(!(rownames(expr_symbol) %in% unlist(c_sets)))]
prots_not_in_sets=rownames(prot_symbol)[which(!(rownames(prot_symbol) %in% unlist(c_sets)))]

sim_prot.prop_x_eff=vector("list")
for(i in 1:length(eff)){
  groups=c(rep("Healthy",50),rep("Case",51))#Set groups
  set.seed(2026)#set seed
  sim_prot.prop_x_eff[[i]]=scenario(list(expr_symbol,prot_symbol),groups,c_sets,n_iter,n_cores,eff[i],0.5,0.25,0)
}

for(i in 1:length(eff)){
  c_sim=sim_prot.prop_x_eff[[i]]
  col_grad=grad(100)[as.numeric(cut(set_info[match(names(c_sets),rownames(set_info)),"size"],100))]#In this case, use set size
  png(paste0("Figures/STsp/prot_prop/",i,"_eff",c(0,0.15,0.4)[i],".png"), units="in", width=11, height=8.5, res=300)
  par(mfrow=c(2,1),mar=c(1.5,4,0,0)+0.5,xpd=t)
  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("DPM_pv")]}),
          ylab="DPM p-value",
          col=col_grad,names=round(mid_cut_set_info$prot.prop,2),ylim=c(0,1))
  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("GSEMA_pv")]}),
          ylab="GSEMA p-value",
          col=col_grad,names=round(mid_cut_set_info$prot.prop,2),ylim=c(0,1))
  par(opar)
  dev.off()
}


names(sim_prot.prop_x_eff)=paste0("eff=",round(eff,2))
save(sim_prot.prop_x_eff,file="Ready/STsp_sim_prot.prop_x_eff.RData")

#TESTING CORRELATION IN SET ----
cut_set_info=data.frame(
  Size=cut(set_info$size,quantile(set_info$size,c(0,0.6,0.85,0.95,1)),dig.lab=2,ordered_result = T),
  prot.prop=cut(set_info$prot.prop,quantile(set_info$prot.prop,c(0,0.65,0.85,0.99,1)),dig.lab=2,ordered_result = T),
  Expr.Cor=cut(set_info$Expr.Cor,quantile(set_info$Expr.Cor,seq(0,1,,11)),dig.lab=2,ordered_result = T),
  Expr.Sd=cut(set_info$Expr.Sd,quantile(set_info$Expr.Sd,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  row.names = rownames(set_info)
)
arr_cont=table(cut_set_info)
apply(arr_cont[2:3,2:3,,2:3], 3, sum)
sum(apply(arr_cont[2:3,2:3,,2:3], 3, sum))#41 sets


mid_levels=t(sapply(cut_set_info[,-3], levels))[,2:3]
mid_cut_set_info=cut_set_info[cut_set_info$Size %in% mid_levels[1,] &
                                cut_set_info$prot.prop %in% mid_levels[2,] &
                                cut_set_info$Expr.Sd %in% mid_levels[3,],]

mid_cut_set_info$Expr.Cor=set_info[row.names(mid_cut_set_info),3]
mid_cut_set_info=mid_cut_set_info[order(mid_cut_set_info$Expr.Cor),]

c_sets=msigdb_C2[rownames(mid_cut_set_info)]

writeLines(
  sapply(lapply(names(c_sets), function(x){c(x,"NA",c_sets[[x]])}), paste, collapse = "\t"),
  con = "Ready/simulation_sets.gmt"
)

genes_not_in_sets=rownames(expr_symbol)[which(!(rownames(expr_symbol) %in% unlist(c_sets)))]
prots_not_in_sets=rownames(prot_symbol)[which(!(rownames(prot_symbol) %in% unlist(c_sets)))]

sim_expr_cor_x_eff=vector("list")
for(i in 1:length(eff)){
  groups=c(rep("Healthy",50),rep("Case",51))#Set groups
  set.seed(2026)#set seed
  c_sim=sim_expr_cor_x_eff[[i]]=scenario(list(expr_symbol,prot_symbol),groups,c_sets,n_iter,n_cores,eff[i],0.5,0.25,0)

  grad=colorRampPalette(c("lightblue","darkblue"))
  col_grad_1=grad(100)[as.numeric(cut(set_info[match(names(c_sets),rownames(set_info)),"size"],100))]
  col_grad_2=grad(100)[as.numeric(cut(set_info[match(names(c_sets),rownames(set_info)),"prot.prop"],100))]


  png(paste0("Figures/STsp/expr_cor/",i,"_eff",c(0,0.15,0.4)[i],".png"), units="in", width=11, height=8.5, res=300)
  par(mfrow=c(2,1),mar=c(5,5,0,5.5)+0.5,xpd=t)
  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("DPM_pv")]}),
          xlab="Pathway inner correlation",ylab="DPM p-value",
          col=col_grad_1,names=round(mid_cut_set_info$Expr.Cor,2),ylim=c(0,1))

  legend("topright",inset=c(-0.13,0),legend=c("size",round(range(set_info[match(names(c_sets),rownames(set_info)),"size"]),2),paste0("eff=",c(0,0.15,0.4)[i])),
         pch=19,col=c("white","lightblue","darkblue","white"),bty="n")

  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("GSEMA_pv")]}),
          xlab="Pathway inner correlation",ylab="GSEMA p-value",
          col=col_grad_2,names=round(mid_cut_set_info$Expr.Cor,2),ylim=c(0,1))

  legend("topright",inset=c(-0.13,0),legend=c("prot.prop",round(range(set_info[match(names(c_sets),rownames(set_info)),"prot.prop"]),2),paste0("eff=",c(0,0.15,0.4)[i])),
         pch=19,col=c("white","lightblue","darkblue","white"),bty="n")


  par(opar)
  dev.off()
}

names(sim_expr_cor_x_eff)=paste0("eff=",round(eff,2))
save(sim_expr_cor_x_eff,file="Ready/STsp_sim_expr_cor_x_eff.RData")

#TESTING DEVIATION IN SET ----
cut_set_info=data.frame(
  Size=cut(set_info$size,quantile(set_info$size,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  prot.prop=cut(set_info$prot.prop,quantile(set_info$prot.prop,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  Expr.Cor=cut(set_info$Expr.Cor,quantile(set_info$Expr.Cor,c(0,0.35,0.5,0.65,1)),dig.lab=2,ordered_result = T),
  Expr.Sd=cut(set_info$Expr.Sd,quantile(set_info$Expr.Sd,seq(0,1,,11)),dig.lab=2,ordered_result = T),
  row.names = rownames(set_info)
)
arr_cont=table(cut_set_info)
apply(arr_cont[2:3,2:3,2:3,], 4, sum)
sum(apply(arr_cont[2:3,2:3,2:3,], 4, sum))


mid_levels=t(sapply(cut_set_info[,-4], levels))[,2:3]
mid_cut_set_info=cut_set_info[cut_set_info$Size %in% mid_levels[1,] &
                                cut_set_info$prot.prop %in% mid_levels[2,] &
                                cut_set_info$Expr.Cor %in% mid_levels[3,],]

mid_cut_set_info$Expr.Sd=set_info[row.names(mid_cut_set_info),4]
mid_cut_set_info=mid_cut_set_info[order(mid_cut_set_info$Expr.Sd),]

c_sets=msigdb_C2[rownames(mid_cut_set_info)]

writeLines(
  sapply(lapply(names(c_sets), function(x){c(x,"NA",c_sets[[x]])}), paste, collapse = "\t"),
  con = "Ready/simulation_sets.gmt"
)

genes_not_in_sets=rownames(expr_symbol)[which(!(rownames(expr_symbol) %in% unlist(c_sets)))]
prots_not_in_sets=rownames(prot_symbol)[which(!(rownames(prot_symbol) %in% unlist(c_sets)))]

sim_expr_sd_x_eff=vector("list")
for(i in 1:length(eff)){
  groups=c(rep("Healthy",50),rep("Case",51))#Set groups
  set.seed(2026)#set seed
  c_sim=sim_expr_sd_x_eff[[i]]=scenario(list(expr_symbol,prot_symbol),groups,c_sets,n_iter,n_cores,eff[i],0.5,0.25,0)

  grad=colorRampPalette(c("lightblue","darkblue"))
  col_grad_1=grad(100)[as.numeric(cut(set_info[match(names(c_sets),rownames(set_info)),"size"],100))]
  col_grad_2=grad(100)[as.numeric(cut(set_info[match(names(c_sets),rownames(set_info)),"prot.prop"],100))]

  png(paste0("Figures/STsp/expr_sd/",i,"_eff",c(0,0.15,0.4)[i],".png"), units="in", width=11, height=8.5, res=300)
  par(mfrow=c(2,1),mar=c(4,4,0,5.5)+0.5,xpd=t)
  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("DPM_pv")]}),
          xlab="Pathway inner standard deviation",ylab="DPM p-value",
          col=col_grad_1,names=round(mid_cut_set_info$Expr.Sd,2),ylim=c(0,1))

  legend("topright",inset=c(-0.13,0),legend=c("size",round(range(set_info[match(names(c_sets),rownames(set_info)),"size"]),2),paste0("eff=",c(0,0.15,0.4)[i])),
         pch=19,col=c("white","lightblue","darkblue","white"),bty="n")

  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("GSEMA_pv")]}),
          xlab="Pathway inner standard deviation",ylab="GSEMA p-value",
          col=col_grad_2,names=round(mid_cut_set_info$Expr.Sd,2),ylim=c(0,1))

  legend("topright",inset=c(-0.13,0),legend=c("prot.prop",round(range(set_info[match(names(c_sets),rownames(set_info)),"prot.prop"]),2),paste0("eff=",c(0,0.15,0.4)[i])),#In this case we mix both
         pch=19,col=c("white","lightblue","darkblue","white"),bty="n")
  par(opar)
  dev.off()
}
names(sim_expr_sd_x_eff)=paste0("eff=",round(eff,2))
save(sim_expr_sd_x_eff,file="Ready/STsp_sim_expr_sd_x_eff.RData")

#SELECTING SETS FOR NEXT SCRIPT ----
#Now we can define the characteristics for the sets used in the script Simulation_Testing_distribution_parameters.R
#DPM starts to not have false positives from size 160 onward, the same goes for proportion of proteins at 0.75% (as we are controling the both parameters, we can take values a bit lower in both), for the other parameters, in order, common values are (retrived from summary) 0.21,0.19,0.35,0.21
#Subtract said values from the set information matrix
temp_set_info=set_info-matrix(c(150,0.7,0.21,0.19,0.36,0.21),nrow=nrow(set_info),ncol=6,T)
#We are interested in positive values for the first 2 variables
temp_set_info=temp_set_info[temp_set_info$size>=0 & temp_set_info$prot.prop>=0,]
#The rest of the variables can be taken at absolute value:
temp_set_info[,3:6]=abs(temp_set_info[,3:6])
#Obtain the first deciles of each variables, and see if some sets fall in the parameters determined before
temp_set_info[]=lapply(temp_set_info,function(x){
  x=cut(x,breaks=quantile(x,seq(0,1,,11)),include.lowest=T,ordered_results=T)})
temp_set_info[temp_set_info$size %in% levels(temp_set_info$size)[1:5] &
                temp_set_info$prot.prop %in% levels(temp_set_info$prot.prop)[2:10] & #Sizes too big can produce a lack of statistical power, but more proportion of proteins is always desirable
                temp_set_info$Expr.Cor %in% levels(temp_set_info$Expr.Cor)[1:5] &
                temp_set_info$Prot.Cor %in% levels(temp_set_info$Prot.Cor)[1:5] &
                temp_set_info$Expr.Sd %in% levels(temp_set_info$Expr.Sd)[1:5] &
                temp_set_info$Prot.Sd %in% levels(temp_set_info$Prot.Sd)[1:5],]

#Now we want a set that has an small size, for the rest of the parameters we can preserve the criteria (except fo proportion of genes present, which we will increase to ensure it does not greatly affect false positives, as we want to isolate the effect from size)

temp_set_info=set_info-matrix(c(0,0.75,0.21,0.19,0.36,0.21),nrow=nrow(set_info),ncol=6,T)
#We are interested in positive values for the first 2 variables
temp_set_info=temp_set_info[temp_set_info$prot.prop>=0,]
#The rest of the variables can be taken at absolute value:
temp_set_info[,3:6]=abs(temp_set_info[,3:6])
#Obtain the first deciles of each variables, and see if some sets fall in the parameters determined before
temp_set_info[]=lapply(temp_set_info,function(x){
  x=cut(x,breaks=quantile(x,seq(0,1,,11)),include.lowest=T,ordered_results=T)})
temp_set_info[temp_set_info$size %in% levels(temp_set_info$size)[1:4] &
                temp_set_info$prot.prop %in% levels(temp_set_info$prot.prop)[1:5] &
                temp_set_info$Expr.Cor %in% levels(temp_set_info$Expr.Cor)[1:5] &
                temp_set_info$Prot.Cor %in% levels(temp_set_info$Prot.Cor)[1:5] &
                temp_set_info$Expr.Sd %in% levels(temp_set_info$Expr.Sd)[1:5] &
                temp_set_info$Prot.Sd %in% levels(temp_set_info$Prot.Sd)[1:5],]



#Now the same for proportion of genes

temp_set_info=set_info-matrix(c(160,0,0.21,0.19,0.36,0.21),nrow=nrow(set_info),ncol=6,T)
#We are interested in positive values for the first 2 variables
temp_set_info=temp_set_info[temp_set_info$size>=0,]
#The rest of the variables can be taken at absolute value:
temp_set_info[,3:6]=abs(temp_set_info[,3:6])
#Obtain the first deciles of each variables, and see if some sets fall in the parameters determined before
temp_set_info[]=lapply(temp_set_info,function(x){
  x=cut(x,breaks=quantile(x,seq(0,1,,11)),include.lowest=T,ordered_results=T)})
temp_set_info[temp_set_info$size %in% levels(temp_set_info$size)[1:5] &
                temp_set_info$prot.prop %in% levels(temp_set_info$prot.prop)[1:3] &
                temp_set_info$Expr.Cor %in% levels(temp_set_info$Expr.Cor)[1:5] &
                temp_set_info$Prot.Cor %in% levels(temp_set_info$Prot.Cor)[1:5] &
                temp_set_info$Expr.Sd %in% levels(temp_set_info$Expr.Sd)[1:5] &
                temp_set_info$Prot.Sd %in% levels(temp_set_info$Prot.Sd)[1:5],]

#Lets check if any of them has some strange behaviour

c_sets=msigdb_C2[c("DEBIASI_APOPTOSIS_BY_REOVIRUS_INFECTION_DN","WELCSH_BRCA1_TARGETS_UP",
                   "BIOCARTA_PITX2_PATHWAY","REACTOME_METABOLISM_OF_COFACTORS","REACTOME_SIGNAL_TRANSDUCTION_BY_L1",
                   "FUJII_YBX1_TARGETS_DN","HORIUCHI_WTAP_TARGETS_DN","PASQUALUCCI_LYMPHOMA_BY_GC_STAGE_DN","VERHAAK_GLIOBLASTOMA_CLASSICAL")]



writeLines(
  sapply(lapply(names(c_sets), function(x){c(x,"NA",c_sets[[x]])}), paste, collapse = "\t"),
  con = "Ready/simulation_sets.gmt"
)

genes_not_in_sets=rownames(expr_symbol)[which(!(rownames(expr_symbol) %in% unlist(c_sets)))]
prots_not_in_sets=rownames(prot_symbol)[which(!(rownames(prot_symbol) %in% unlist(c_sets)))]

sim_test_sets=vector("list",length = 3)
eff=qnorm(0.5+c(0,0.15,0.4)/2)
for(i in 1:3){
  groups=c(rep("Healthy",50),rep("Case",51))#Set groups
  set.seed(2026)#set seed
  c_sim=sim_test_sets=scenario(list(expr_symbol,prot_symbol),groups,c_sets,n_iter,n_cores,eff[i],0.5,0.25,0)

  png(paste0("Figures/STsp/test_sets/",i,"_eff",c(0,0.15,0.4)[i],".png"), units="in", width=11, height=8.5, res=300)
  par(mfrow=c(2,1),mar=c(2,2,1,0)+0.5)
  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("DPM_pv")]}),
          main=paste0("DPM p-value distribution by std. deviation in set for effect=",round(eff[i],2)),
          col=c("steelblue","steelblue","goldenrod","goldenrod","goldenrod","indianred2","indianred2","indianred2","indianred2"),
          names=c("C11","C12","C21","C22","C23","C31","C32","C33","C34"),ylim=c(0,1))

  boxplot(sapply(c_sim[match(names(c_sets),names(c_sim))], function(x){x[,c("GSEMA_pv")]}),
          main=paste0("GSEMA p-value distribution by std. deviation in set for effect=",round(eff[i],2)),
          col=c("steelblue","steelblue","goldenrod","goldenrod","goldenrod","indianred2","indianred2","indianred2","indianred2")
          ,names=c("C11","C12","C21","C22","C23","C31","C32","C33","C34"),ylim=c(0,1))
  dev.off()
}

#We can take the first one of each category
c_sets=c_sets[c(1,3,6)]#This sets are c("DEBIASI_APOPTOSIS_BY_REOVIRUS_INFECTION_DN","BIOCARTA_PITX2_PATHWAY","FUJII_YBX1_TARGETS_DN")
save("c_sets",file="Ready/simulation_dist.RData")

