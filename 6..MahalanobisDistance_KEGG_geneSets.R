################################### Array data ###################################
rm(list=ls())
library(pracma)
setwd("/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/")
source("SCRIPTS/FinalScripts/supplementary/Functions.R")

data=read.csv("DATA/Voineagu2011_Array.csv") # read in data
kegg=read.csv("DATA/c2.cp.kegg.v6.2.symbols.csv")

# select gene sets with at least half of the genes present in the dataset
l=rep(0, ncol(kegg))
for (j in c(1:length(l))) l[j]=length(intersect(data$Symbol, kegg[,j]))/(nrow(kegg)-length(which(kegg[,j]%in%"")))
kegg=kegg[, which(l>0.5)]
symbols=data$Symbol
data=data[,-c(1,2)] # remove gene annotation columns

results_kegg_MD_Array=run_MD_geneSets(expData=data,symbols=symbols, kegg=kegg,num_trials=100, transposed=TRUE, asd_label="A_", con_label="C_")
names(results_kegg_MD_Array)=colnames(kegg)
results_kegg_MD_Array[which(results_kegg_MD_Array< 0.01)]
save(results_kegg_MD_Array, file="RESULTS//FINAL/Mahalanobis_Kegg/results_kegg_MD_Array.rda")


################################### Replication in RNA-seq data ###################################
rm(list=ls())

setwd("/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/")
source("SCRIPTS/FinalScripts/supplementary/Functions.R")
data=read.csv("DATA/Parikshak2016_RnaSeq.csv") # read in data
kegg=read.csv("DATA/c2.cp.kegg.v6.2.symbols.csv")
load("RESULTS//FINAL/Mahalanobis_Kegg/results_kegg_MD_Array.rda")
kegg=kegg[,which(colnames(kegg)%in%names(results_kegg_MD_Array[which(results_kegg_MD_Array < 0.01)])) ]


hgnc=read.table("/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/ANNOTATION/Ensembl_to_other/HGNC_EnsemblGeneID.txt", sep="\t", header=TRUE)
symbols=hgnc$Approved.Symbol[match(data[,1], hgnc$EnsemblGeneID)]
data=data[,-1] # remove gene annotation columns

results_kegg_MD_RNAseq=run_MD_geneSets(expData=data,symbols=symbols, kegg=kegg,num_trials=100, transposed=TRUE, asd_label="ASD_", con_label="CTL_")
save(results_kegg_MD_RNAseq, file="RESULTS/FINAL/Mahalanobis_Kegg/results_kegg_MD_RNAseq.rda")

