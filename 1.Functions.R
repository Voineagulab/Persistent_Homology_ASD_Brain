#########################################################################################################
SumDeathTime=function(diagram)
#Function that calculates the Sum of Death Times of 0-dimensional componenents of a persistent homology diagram
{
  diagram=as.data.frame(diagram)
  colnames(diagram)=c("Dim", "Birth", "Death")
  index=which(diagram$Dim==0)
  return(sum(diagram[index,3]))
}
#########################################################################################################
EulerChar=function(diagram)
#Function that calculates the Euler Characteristic of a persistent homology diagram
{ 
  diagram=as.data.frame(diagram)
  colnames(diagram)=c("Dim", "Birth", "Death")
  diff=diagram$Death-diagram$Birth
  return(sum(diff[which(diagram$Dim==0)]) -
           sum(diff[which(diagram$Dim==1)]) +
           sum(diff[which(diagram$Dim==2)]))
}
#########################################################################################################
runGroup_pHom=function(ASD, health, num_trials, analysis_type)
{ library(phom)
# This is  a wrapper funtion that compares persistent homology diagrams of two groups of gene expression data.
# It calculates the difference in Sum-of-death times and Euler characteristic between the two diagrams, and assesses statistical significance by random permutation of group labels
  # Parameters:
  # ASD:expression data from ASD samples; health: expression data from control samples.
  # num_trials: number of random samplings 
  # analysis_type: either "bySample" or "byGene" - specifies whether the analysis is done across samples or across genes.
  ###STEP1
  # Calculate Distance Matrices for ASD and control data 
  print ("Working on : pHom analysis of ASD and Controls")
  d.ASD.real=as.matrix(1-cor(ASD, method="p"))
  d.health.real=as.matrix(1-cor(health, method="p"))
  # Calculate Persistent Diagrams
  # For analysis "by Sample" Vietori-Ripps complexes are built (mode = "vr").
  # For analysis "byGene", where the size of the data is very large, the lazy-witness option is used (mode = "lw").
  # Since higher dimensional components have not been observed for the analysis "bySample" of our datasets, the dimension is set to 2 in this case to reduce computational time.
  # if (analysis_type=="bySample")
  # {
  ASD.persistance.real=pHom(d.ASD.real,dimension=2,max_filtration_value = max(d.ASD.real),mode = "vr",metric = "distance_matrix")
  health.persistance.real=pHom(d.health.real,dimension=2,max_filtration_value = max(d.health.real),mode = "vr",metric = "distance_matrix")
  # }
  # if (analysis_type=="byGene")
  # {
  # ASD.persistance.real=pHom(d.ASD.real,dimension=3,max_filtration_value = max(d.ASD.real),mode = "lw",metric = "distance_matrix" , landmark_set_size =20)
  # health.persistance.real=pHom(d.health.real,dimension=3,max_filtration_value = max(d.health.real),mode = "lw",metric = "distance_matrix", landmark_set_size =20)
  # }
  # Calculate the difference in Sum of Death Time and Euler Characteristic between ASD and Control diagrams
  real.diff.sdt=SumDeathTime(ASD.persistance.real)-SumDeathTime(health.persistance.real)
  real.diff.euler=EulerChar(ASD.persistance.real)-EulerChar(health.persistance.real)
  ### STEP2 random permutation of group assignment of samples
  print ("Working on : Random Permutations")
  #placeholder vectors which will contain the difference in Sum of Death Time and Euler Characteristic between ASD and Control diagrams at each random permutation
  trials.diff.sdt=rep(0,num_trials)
  trials.diff.euler=rep(0,num_trials)
  
  #random permutations of group labels
  expData=cbind(ASD,health)
  index=1:(ncol(ASD) + ncol(health))
  
  for (i in 1:num_trials) 
  {
    index=sample(index)
    ASD.trial=expData[,index[1:ncol(ASD)] ]
    health.trial=expData[,index[(ncol(ASD)+1):ncol(expData)] ]
    ## Calculate distance matrix at each random sampling
    d.ASD.trial=as.matrix(1-cor(ASD.trial, method="p"))
    d.health.trial=as.matrix(1-cor(health.trial, method="p"))
    ## Calculate persistent diagrams at each random sampling
   # if (analysis_type=="bySample")
   #  {
    ASD.persistance.trial=pHom(d.ASD.trial,dimension=2,max_filtration_value = max(d.ASD.trial),mode = "vr",metric = "distance_matrix")
    health.persistance.trial=pHom(d.health.trial,dimension=2,max_filtration_value = max(d.health.trial),mode = "vr",metric = "distance_matrix")
    # }
    # if (analysis_type=="byGene")
    # {
    # ASD.persistance.trial=pHom(d.ASD.trial,dimension=3,max_filtration_value = max(d.ASD.trial),mode = "lw",metric = "distance_matrix" , landmark_set_size =20)
    # health.persistance.trial=pHom(d.health.trial,dimension=3,max_filtration_value = max(d.health.trial),mode = "lw",metric = "distance_matrix", landmark_set_size =20)
    # }
    ##retain the difference in Sum of death Time and Euler Characteristic between ASD and controls at each random sampling
    trials.diff.sdt[i]=SumDeathTime(ASD.persistance.trial)-SumDeathTime(health.persistance.trial)
    trials.diff.euler[i]=EulerChar(ASD.persistance.trial)-EulerChar(health.persistance.trial)
    print(i)
  }
  
###STEP3 compare real and trial diff.sdt and diff.euler values 
  print ("Working on : Calculating p-values")
  p.sdt=length(which(trials.diff.sdt > real.diff.sdt))/num_trials
  p.euler=length(which(trials.diff.euler > real.diff.euler))/num_trials
  
###OUTPUT RESULTS
  results=list(ASD.persistance.real=ASD.persistance.real,
               health.persistance.real=health.persistance.real, 
               real.diff.sdt= real.diff.sdt,
               real.diff.euler= real.diff.euler,
               trials.diff.sdt=trials.diff.sdt,
               trials.diff.euler=trials.diff.euler,
               p.sdt=p.sdt,
               p.euler=p.euler)
  return(results)
  print("DONE")
}
#########################################################################################################
# This function carries out the pHom analysis as above, using genes in a specific gene set (G)
# The difference in SDT and Euler Characteristic between ASD and controls is compared to randomly sampled gene sets of of the same size as G.
# Parameters:
# ASD: gene expression data for genes in gene set G, in ASD samples
# health: gene expression data for genes in gene set G, in Control samples
# data: gene xpression data for all genes and all samples (used for random sampling of genes)
# num_trials: number of random samplings
# asd_label and control_label: group-specific prefix of sample names.
runGroup_pHom_geneSets=function(ASD, health, data, num_trials, asd_label, con_label)
{ library(phom)
  print ("Working on : geneSet")
  d.ASD.real=as.matrix(1-cor(ASD, method="p"))
  d.health.real=as.matrix(1-cor(health, method="p"))
  
  ASD.persistance.real=pHom(d.ASD.real,dimension=2,max_filtration_value = max(d.ASD.real),mode = "vr",metric = "distance_matrix")
  health.persistance.real=pHom(d.health.real,dimension=2,max_filtration_value = max(d.health.real),mode = "vr",metric = "distance_matrix")
  
  real.diff.sdt=SumDeathTime(ASD.persistance.real)-SumDeathTime(health.persistance.real)
  real.diff.euler=EulerChar(ASD.persistance.real)-EulerChar(health.persistance.real)
  ### STEP2 #random sampling of gene sets of the same size as G.
  #placeholder vectors 
  trials.diff.sdt=rep(0,num_trials)
  trials.diff.euler=rep(0,num_trials)
  #random sampling
  for (i in 1:num_trials) 
  {
    s=sample(c(1:nrow(data)) , nrow(ASD))
    ASD.trial=data[s,grep(asd_label, colnames(data))]
    health.trial=data[s,grep(con_label, colnames(data))]
    ## Calculate distance matrix at each random sampling
    d.ASD.trial=as.matrix(1-cor(ASD.trial, method="p"))
    d.health.trial=as.matrix(1-cor(health.trial, method="p"))
    ## Calculate persistent diagrams at each random sampling
    # if (analysis_type=="bySample")
    #  {
    ASD.persistance.trial=pHom(d.ASD.trial,dimension=2,max_filtration_value = max(d.ASD.trial),mode = "vr",metric = "distance_matrix")
    health.persistance.trial=pHom(d.health.trial,dimension=2,max_filtration_value = max(d.health.trial),mode = "vr",metric = "distance_matrix")
    # }
    # if (analysis_type=="byGene")
    # {
    # ASD.persistance.trial=pHom(d.ASD.trial,dimension=3,max_filtration_value = max(d.ASD.trial),mode = "lw",metric = "distance_matrix" , landmark_set_size =20)
    # health.persistance.trial=pHom(d.health.trial,dimension=3,max_filtration_value = max(d.health.trial),mode = "lw",metric = "distance_matrix", landmark_set_size =20)
    # }
    ##retain the difference in Sum of death Time and Euler Characteristic between ASD and controls at each random sampling
    trials.diff.sdt[i]=SumDeathTime(ASD.persistance.trial)-SumDeathTime(health.persistance.trial)
    trials.diff.euler[i]=EulerChar(ASD.persistance.trial)-EulerChar(health.persistance.trial)
    print(i)
  }
  
  ###STEP3 compare real and trial diff.sdt and diff.euler values 
  print ("Working on : Calculating p-values")
  p.sdt=length(which(trials.diff.sdt > real.diff.sdt))/num_trials
  p.euler=length(which(trials.diff.euler > real.diff.euler))/num_trials
  
  ###OUTPUT RESULTS
  results=list(p.sdt=p.sdt,p.euler=p.euler)
  return(results)
  print("DONE")
}

#########################################################################################################
# This function carries out a sample-level Mahalanobis distance analysis
# Parameters:
# ASD: gene expression data for ASD samples
# ASD: gene expression data for Control samples
# expData: gene expression data for all samples
# num_trials: number of random samplings
# asd_label and control_label: group-specific prefix of sample names.
run_MD=function( ASD, CON,  expData, num_trials, asd_label, con_label)
{
library(pracma)
##### Calculate Sum of Squared Mahalanobis Distances (SSMD) for ASD samples relative to CONTROLS
# generalized inverse of the covariance matrix of controls
cov.inv=pinv(cov(CON))
# Squared Mahalanobis Distance
SMD=mahalanobis(ASD, apply(CON,2,mean), cov.inv, inverted=TRUE)
# Sum of Squared Mahalanobis Distance
SSMD=sum(SMD)
##### Random permutation
trials.SSMD = rep(0,num.trials)
index=1:(nrow(ASD)+nrow(CON))

for (i in 1:num.trials) 
{ print(i)
  index=sample(index)
  ASD.test=t(expData[,index[ASD_samples]])
  CON.test=t(expData[,index[CON_samples]])
  
  cov.inv.test=pinv(cov(CON.test))
  SMD.test=mahalanobis(ASD.test, apply(CON.test,2,mean), cov.inv.test, inverted=TRUE)
  trials.SSMD[i]=sum(SMD.test)
}
results=list(SSMD=SSMD, trials.SSMD=trials.SSMD)
return(results)
}
#########################################################################################################
# This function carries out a Mahalanobis distance analysis as above, using genes in a specific gene set (G)
# The SSMD value for ASD samples is compared to randomly sampled gene sets of of the same size as G.
# Parameters:
# expData: gene expression data; Note that the selection of samples in ASD and control groups, and genes in gene set G is done within the function
# symbols: gene symbols for the gene expression dataset
# kegg: table of KEGG gene sets, with one gene set per column.
# num_trials: number of random samplings
# asd_label and control_label: group-specific prefix of sample names.
# 
run_MD_geneSets=function(expData,symbols, kegg,num_trials, transposed, asd_label, con_label)
{
library(pracma)
results_kegg_md=rep(NA, ncol(kegg))
######### MD on each pathway
for (j in c(1:ncol(kegg)))
{
  print(paste("----", j, colnames(kegg)[j]))
  k=which(symbols%in%kegg[,j])
  ASD=t(expData[k,grep(asd_label, colnames(expData))]) # select ASD samples
  CON=t(expData[k,grep(con_label, colnames(expData))]) # select control samples

  cov.inv=pinv(cov(CON))
  # Squared Mahalanobis Distance
  SMD=mahalanobis(ASD, apply(CON,2,mean), cov.inv, inverted=TRUE)
  # Sum of Squared Mahalanobis Distance
  SSMD=sum(SMD)
  #random sampling of gene sets of the same size as G
  trials.SSMD=rep(NA,num_trials)
  for (i in 1:num_trials) 
  { #print(i)
    s=sample(c(1:nrow(expData)) , length(k))
    ASD.trial=t(expData[s,grep(asd_label, colnames(expData))])
    CON.trial=t(expData[s,grep(con_label, colnames(expData))])
    
    cov.inv.trial=pinv(cov(CON.trial))
    # Squared Mahalanobis Distance
    SMD.trial=mahalanobis(ASD.trial, apply(CON.trial,2,mean), cov.inv.trial, inverted=TRUE)
    # Sum of Squared Mahalanobis Distance
    trials.SSMD[i]=sum(SMD.trial)
  }
  results_kegg_md[j]=length(which(trials.SSMD > SSMD))/length(trials.SSMD)
}
names(results_kegg_md)=colnames(kegg)
return(results_kegg_md)
}
