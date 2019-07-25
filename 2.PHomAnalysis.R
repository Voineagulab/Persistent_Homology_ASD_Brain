
################################### Array data ###################################
rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/")
source("SCRIPTS/FinalScripts/supplementary/Functions.R")
data=read.csv("DATA/Voineagu2011_Array.csv") # read in data
data=data[,-c(1,2)] # remove gene annotation columns
ASD=data[,grep("A_", colnames(data))] # select ASD samples
health=data[,grep("C_", colnames(data))] # select control samples
######### pHom analysis By Samples
results_Array_bySample=runGroup_pHom(ASD=ASD, health=health, analysis_type = "bySample", num_trials=100000)
save(results_Array_bySample, file = "RESULTS/FINAL/results_Array_bySample_100000.rda" )
results_Array_bySample=runGroup_pHom(ASD=ASD, health=health, analysis_type = "bySample", num_trials=1000)
save(results_Array_bySample, file = "RESULTS/FINAL/results_Array_bySample_1000.rda" )
# ######### pHom analysis By Genes
results_Array_byGene=runGroup_pHom(ASD=t(ASD), health=t(health), analysis_type = "byGene", num_trials=1000)
save(results_Array_byGene, file = "RESULTS/FINAL/results_Array_byGene.rda_1000" )


################################### RNA-seq data ###################################
rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/")
source("SCRIPTS/FinalScripts/supplementary/Functions.R")
data=read.csv("DATA/Parikshak2016_RnaSeq.csv") # read in data
data=data[,-1] # remove gene annotation columns
ASD=data[,grep("ASD_", colnames(data))] # select ASD samples
health=data[,grep("CTL_", colnames(data))] # select control samples
################################### pHom analysis By Samples
results_RNAseq_bySample=runGroup_pHom(ASD=ASD, health=health,  analysis_type = "bySample", num_trials=1000)
save(results_RNAseq_bySample, file = "RESULTS/FINAL/results_RNAseq_bySample_1000.rda" )

# ################################### pHom analysis By Genes
results_RNAseq_byGene=runGroup_pHom(ASD=t(ASD), health=t(health),  file_name="RNAseq_byGene", analysis_type = "byGene", num_trials=1000)
save(results_RNAseq_byGene, file = "RESULTS/FINAL/results_RNAseq_byGene_1000.rda" )



