################################### Array data ###################################
rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/")
source('/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/SCRIPTS/FinalScripts/supplementary/Functions.R')

##### Load in and format data
data=read.csv("DATA/Voineagu2011_Array.csv")
expData=data[,-c(1,2)] # remove gene annotation columns
ASD=t(expData[,grep("A_", colnames(expData))])
CON=t(expData[,grep("C_", colnames(expData))])
##### Run MD
results.Mahalanobis.array=run_MD( ASD=ASD, CON=CON, num_trials=1000, expData=expData, asd_label="A_", con_label="C_")
##### Save
save(results.Mahalanobis.array, file="RESULTS/FINAL/Mahalanobis/results.Mahalanobis.array_1000.rda")

################################### RNA-seq data ###################################
rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/")
source('/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/SCRIPTS/FinalScripts/supplementary/Functions.R')

##### Load in and format data
data=read.csv("DATA/Parikshak2016_RnaSeq.csv") # read in data
data=data[,-1] # remove gene annotation columns
ASD=t(data[,grep("ASD_", colnames(data))]) # select ASD samples
CON=t(data[,grep("CTL_", colnames(data))]) # select control samples

##### Run MD
start=Sys.time()
results.Mahalanobis.array=run_MD( ASD=ASD, CON=CON, num_trials=1, expData=expData, asd_label="A_", con_label="C_")
end=Sys.time()
time=end-start
save(time,file="RESULTS/FINAL/MD.RNA-seq.time.rda")
