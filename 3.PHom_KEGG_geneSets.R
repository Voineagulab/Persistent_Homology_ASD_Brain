################################### Array data ###################################
rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/")
source("SCRIPTS/FinalScripts/supplementary/Functions.R")
data=read.csv("DATA/Voineagu2011_Array.csv") # read in data
kegg=read.csv("DATA/c2.cp.kegg.v6.2.symbols.csv")
l=rep(0, ncol(kegg))
for (j in c(1:length(l))) l[j]=length(intersect(data$Symbol, kegg[,j]))/(nrow(kegg)-length(which(kegg[,j]%in%"")))
kegg=kegg[, which(l>0.5)]
symbols=data$Symbol

data=data[,-c(1,2)] # remove gene annotation columns
ASD=data[,grep("A_", colnames(data))] # select ASD samples
health=data[,grep("C_", colnames(data))] # select control samples

results_kegg_sdt=rep(NA, ncol(kegg))
results_kegg_euler=rep(NA, ncol(kegg))
######### pHom analysis on each pathway
for (j in c(1:ncol(kegg)))
{
print(paste("----", j, colnames(kegg)[j]))
k=which(symbols%in%kegg[,j])
results=runGroup_pHom_geneSets(ASD=ASD[k,], health=health[k,], data,num_trials=100, asd_label="A_", con_label="C_")
results_kegg_sdt[j]=results$p.sdt
results_kegg_euler[j]=results$p.euler
}
names(results_kegg_sdt)=colnames(kegg)
save(results_kegg_sdt, file="RESULTS/FINAL/PHom_Kegg/results_kegg_sdt_Array.rda")
names(results_kegg_euler)=colnames(kegg)
save(results_kegg_euler, file="RESULTS/FINAL/PHom_Kegg/results_kegg_euler_Array.rda")

################################### Replication in RNA-seq data ###################################
rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/")
source("SCRIPTS/FinalScripts/supplementary/Functions.R")
data=read.csv("DATA/Parikshak2016_RnaSeq.csv") # read in data
kegg=read.csv("DATA/c2.cp.kegg.v6.2.symbols.csv")
load(file="RESULTS/FINAL/PHom_Kegg/results_kegg_sdt_Array.rda")
load(file="RESULTS/FINAL/PHom_Kegg/results_kegg_euler_Array.rda")

sig=which((results_kegg_euler < 0.01) & (results_kegg_sdt < 0.01))
kegg=kegg[,c(1,which(colnames(kegg)%in% names(results_kegg_sdt)[sig]))]

hgnc=read.table("/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/ANNOTATION/Ensembl_to_other/HGNC_EnsemblGeneID.txt", sep="\t", header=TRUE)
symbols=hgnc$Approved.Symbol[match(data[,1], hgnc$EnsemblGeneID)]
data=data[,-1] # remove gene annotation columns

ASD=data[,grep("ASD_", colnames(data))] # select ASD samples
health=data[,grep("CTL_", colnames(data))] # select control samples

results_kegg_sdt=rep(NA, ncol(kegg))
results_kegg_euler=rep(NA, ncol(kegg))

######### pHom analysis on each pathway
for (j in c(2:ncol(kegg)))
{
  print(paste("----", j, colnames(kegg)[j]))
  k=which(symbols%in%kegg[,j])
  results=runGroup_pHom_geneSets(ASD=ASD[k,], health=health[k,], data,num_trials=100, asd_label="ASD_", con_label="CTL_")
  results_kegg_sdt[j]=results$p.sdt
  results_kegg_euler[j]=results$p.euler
}
save(results_kegg_sdt, file="RESULTS//FINAL/PHom_Kegg/results_kegg_sdt_RNAseq.rda")
save(results_kegg_euler, file="RESULTS//FINAL/PHom_Kegg/results_kegg_euler_RNAseq.rda")
