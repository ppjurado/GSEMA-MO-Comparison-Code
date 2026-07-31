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
#Set simulation parameters
n_iter=100
n_cores=10


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
    dpm_pval_matrix, gmt = "Ready/simulation_dist.gmt",significant = 1,correction_method = "none", merge_method = "DPM",
    scores_direction = dpm_dir_matrix, constraints_vector = dpm_constraints_vector, geneset_filter = c(5, 10000))
  return(as.data.frame(dpm_results[,c(1,3)]))#Return path name and nominal p-value (parameter correction_method = "none")
}

iteration=function(data,group,sets){#Join the results
  dpm_res=iteration_DPM(data,group,sets)
  GSEMA_res=iteration_GSEMA(data,group,sets)
  it_res=cbind(dpm_res,GSEMA_res)
  colnames(it_res)=c("Path","DPM_pv","GSEMA_eff","GSEMA_eff_var","GSEMA_pv")
  return(it_res)
}

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

#Load the genesets determined in Simulation_Testing_set_parameters.R
load("Ready/simulation_dist.RData")

names(c_sets)=c("S1","S2","S3")

writeLines(
  sapply(lapply(names(c_sets), function(x){c(x,"NA",c_sets[[x]])}), paste, collapse = "\t"),
  con = "Ready/simulation_dist.gmt"
)

#Load LUAD data (saved in LUAD_Enrichment_MO_Analisys)
load("Ready/LUAD_expr_prot.RData")

#Simulations will be performed only with controls, as they have sufficient samples + the matrix for simulation will be more homogeneus
expr_symbol=exprNormal_symbol[-which(rowSums(exprNormal_symbol)==0),]
prot_symbol=protNormal_symbol
rm(exprTumor_symbol,protTumor_symbol,exprNormal_symbol,protNormal_symbol)
#See which genes and proteins are outside the sets
genes_not_in_sets=rownames(expr_symbol)[which(!(rownames(expr_symbol) %in% unlist(c_sets)))]
prots_not_in_sets=rownames(prot_symbol)[which(!(rownames(prot_symbol) %in% unlist(c_sets)))]
#Simulations will mostly vary in the parameter effect

#Positive effect sizes ----
#Define the values that parameters will be iterated over:
sim_param_comb=expand.grid(seq(0,0.60,by=0.06),#11 effect sizes (as percentage now for naming purposes, later will change to actual normal value)
                           c(0.25,0.5,0.75),#3 proportions of genes affected
                           c(0.25),
                           c(1))

rownames(sim_param_comb)=apply(round(sim_param_comb,2),1,function(x) {sprintf("F%0.2f_Ps%0.2f_Pb%0.2f",x[1],x[2],x[3])})
sim_param_comb[,1]=qnorm(0.5+sim_param_comb[,1]/2)


sim_mmp=apply(sim_param_comb, 1, function(x){
  groups=c(rep("Healthy",50),rep("Case",51))#Set groups
  set.seed(2025)
  scenario(list(expr_symbol,prot_symbol),groups,c_sets,n_iter,n_cores,x[1],x[2],x[3],x[4])
})
#sim_mmp contains a list, each element of the list is a combination of distributional parameters (scenario), inside each, there is another list containing each of the sets, inside each set there are (in matrices) the results of both methods (columns) by each iteration (rows)
#now form a long format matrix with the p-values of each method

#Change name of column sim for ordering in boxplot
name_mapping=data.frame(old=names(sim_mmp),
                        new=apply(round(sim_param_comb,2),1,function(x) {sprintf("Ps%0.2f_F%0.2f_Pb%0.2f",x[2],x[1],x[3])}))

names(sim_mmp)=name_mapping$new

sim_res_DPM=do.call(rbind,#do.call binds the rows into the long format
                    lapply(names(sim_mmp), function(x){#this lapply extracts, for each iteration, a data frame
                      data.frame(#We input in a data frame
                        do.call(rbind,#The combination of rows
                                lapply(sim_mmp[[x]], function(x){#Of each element of a scenario, that is, sim_mmp[[x]] contains for a given set of parameters, a list with elements the gen sets, inside each element, the results of the iterations, and it ads a colum with the name of the set (for ease of plotting)
                        data.frame(pv=x[,1],set=rownames(x)[1])})),
                        sim=x)#Add the name of the scenario
                    }))#The result is a long format matrix of the p-values with their associated pathway and scenario for DPM
#now the same with GSEMA
sim_res_GSEMA=do.call(rbind,
                     lapply(names(sim_mmp), function(x){
                       data.frame(do.call(rbind,lapply(sim_mmp[[x]], function(x){
                         data.frame(pv=x[,4],set=rownames(x)[1])})),sim=x)
                     }))

pdf("Figures/STdp/Jurado-Bascon_Fig3.pdf", width = 13, height = 9,bg="white")
#png(paste0("Figures/STdp/design1.png"), units="in", width=11, height=8.5, res=300)
par(mfrow=c(2,1),mar=c(2.5,4,0,0)+0.5,xpd=T)
boxplot(pv~sim:set,data=sim_res_DPM,
        at=rep(seq(0, by=nrow(sim_param_comb)+1,length.out=length(c_sets)),each=nrow(sim_param_comb))+1:nrow(sim_param_comb),
        las=2,cex.axis=0.7,
        xlab="",ylab="DPM p-value",ylim=c(0,1),
        names=rep(sub(".*?(F[0-9.]+).*","\\1",rownames(sim_param_comb)),3),
        col=rep(c("brown3","darkorange3","goldenrod2"),each=11))
boxplot(pv~sim:set,data=sim_res_GSEMA,
        at=rep(seq(0, by=nrow(sim_param_comb)+1,length.out=length(c_sets)),each=nrow(sim_param_comb))+1:nrow(sim_param_comb),
        las=2,cex.axis=0.7,
        xlab="",ylab="GSEMA p-value",ylim=c(0,1),
        names=rep(sub(".*?(F[0-9.]+).*","\\1",rownames(sim_param_comb)),3),
        col=rep(c("brown3","darkorange3","goldenrod2"),each=11))
legend("topright",inset=c(0,0),legend=c("Ps=0.25","Ps=0.5","Ps=0.75"),pch=19,col=c("brown3","darkorange3","goldenrod2"),bty="n",cex=1)
par(opar)
dev.off()


#DISAGREEMENT IN SETS ----
#p-value merging methods, as pointed by DPM, do not usually take into account the directionality of the p-values, DPM does but marginally, gene-protein pair one by one, it does not take into account if a set is over or infra expressed


sim_param_comb=expand.grid(c(0,0.1,0.25,0.6,0.95),
                           c(0.5),
                           c(0.25),
                           seq(0,1,,5))

rownames(sim_param_comb)=apply(round(sim_param_comb,2),1,function(x){sprintf("Pd%0.2f_F%0.2f_Ps%0.2f_Pb%0.2f",x[4],x[1],x[2],x[3])})
sim_param_comb[,1]=qnorm(0.5+sim_param_comb[,1]/2)

sim_dis=apply(sim_param_comb, 1, function(x){
  groups=c(rep("Healthy",50),rep("Case",51))#Set groups
  set.seed(2026)#set seed
  scenario(list(expr_symbol,prot_symbol),groups,c_sets,n_iter,n_cores,x[1],x[2],x[3],x[4])
})


sim_res_DPM=do.call(rbind,
                    lapply(names(sim_dis), function(x){
                      data.frame(do.call(rbind,lapply(sim_dis[[x]], function(x){
                        data.frame(pv=x[,1],set=rownames(x)[1])})),sim=x)
                    }))

sim_res_GSEMA=do.call(rbind,
                      lapply(names(sim_dis), function(x){
                        data.frame(do.call(rbind,lapply(sim_dis[[x]], function(x){
                          data.frame(pv=x[,4],set=rownames(x)[1])})),sim=x)
                      }))

pdf("Figures/STdp/Jurado-Bascon_Fig4.pdf", width = 13, height = 9,bg="white")
#png(paste0("Figures/STdp/design2.png"), units="in", width=11, height=8.5, res=300)
par(mfrow=c(2,1),mar=c(2.5,4,0,0)+0.5,xpd=T)
boxplot(pv~sim:set,data=sim_res_DPM,
        at=rep(seq(0, by=nrow(sim_param_comb)+1,length.out=length(c_sets)),each=nrow(sim_param_comb))+1:nrow(sim_param_comb),
        las=2,cex.axis=0.75,
        xlab="",ylab="DPM p-value",ylim=c(0,1),
        names=rep(sub(".*?(F[0-9.]+).*","\\1",rownames(sim_param_comb)),3),
        col=rep(rep(RColorBrewer::brewer.pal(5,"PRGn"),each=5),3))
boxplot(pv~sim:set,data=sim_res_GSEMA,
        at=rep(seq(0, by=nrow(sim_param_comb)+1,length.out=length(c_sets)),each=nrow(sim_param_comb))+1:nrow(sim_param_comb),
        las=2,cex.axis=0.75,
        xlab="",ylab="GSEMA p-value",ylim=c(0,1),
        names=rep(sub(".*?(F[0-9.]+).*","\\1",rownames(sim_param_comb)),3),
        col=rep(rep(RColorBrewer::brewer.pal(5,"PRGn"),each=5),3))
legend("topright",inset=c(0,0),legend=c("Pd=0","Pd=0.3","Pd=0.5","Pd=0.3","Pd=1"),pch=19,col=RColorBrewer::brewer.pal(5,"PRGn"),bty="n",cex=0.75)
par(opar)
dev.off()


#save both simulations
save("sim_mmp","sim_dis",file="Ready/STdp.RData")
