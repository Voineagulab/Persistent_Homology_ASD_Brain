setwd("/Volumes/Data1/PROJECTS/PersistentHomology_ASD_Brain/")
rm(list=ls())
################################### Plot Figures ###################################
library(phom)
par(family="serif", cex=1)

###################################Figure3 
pdf("RESULTS/FINAL/Figure_bySample.pdf", height=7.5, width=15)
par(mfrow=c(2,6))
load("RESULTS/FINAL/results_Array_bySample_100000.rda")
results=results_Array_bySample
n_trials=length(results$trials.diff.sdt)
p.sdt=length(which(results$trials.diff.sdt >results$real.diff.sdt))/n_trials
p.euler=length(which(results$trials.diff.euler >results$real.diff.euler))/n_trials

plotPersistenceDiagram(results$ASD.persistance.real, 2,max(results$ASD.persistance.real[,3]), title="Microrray ASD")
plotPersistenceDiagram(results$health.persistance.real, 2, max(results$health.persistance.real[,3]),title="Microrray Controls")

plotPersistenceDiagram(results$ASD.persistance.real, 2,max(results$ASD.persistance.real[which(results$ASD.persistance.real[,1]>0),3]),title="Microrray ASD")
plotPersistenceDiagram(results$health.persistance.real, 2,max(results$health.persistance.real[which(results$health.persistance.real[,1]>0),3]),title="Microrray Controls")

plot(density(results$trials.diff.sdt), sub=round(p.sdt, 5), xlab="Diff.SDT", main= "SDT");abline(v=results$real.diff.sdt, col="red")
plot(density(results$trials.diff.euler), sub=round(p.euler, 5), xlab="Diff.Euler", main= "Euler Characteristic");abline(v=results$real.diff.euler, col="red")

load("RESULTS/FINAL/results_RNAseq_bySample_1000.rda")
results=results_RNAseq_bySample
n_trials=length(results$trials.diff.sdt)
p.sdt=length(which(results$trials.diff.sdt >results$real.diff.sdt))/n_trials
p.euler=length(which(results$trials.diff.euler >results$real.diff.euler))/n_trials

plotPersistenceDiagram(results$ASD.persistance.real, 2,max(results$ASD.persistance.real[,3]), title="RNA-seq ASD")
plotPersistenceDiagram(results$health.persistance.real, 2, max(results$health.persistance.real[,3]),title="RNA-seq Controls")

plotPersistenceDiagram(results$ASD.persistance.real, 2,max(results$ASD.persistance.real[which(results$ASD.persistance.real[,1]>0),3]),title="RNA-seq ASD")
plotPersistenceDiagram(results$health.persistance.real, 2,max(results$health.persistance.real[which(results$health.persistance.real[,1]>0),3]),title="RNA-seq Controls")

plot(density(results$trials.diff.sdt), sub=round(p.sdt, 5), xlab="Diff.SDT", main= "SDT");abline(v=results$real.diff.sdt, col="red")
plot(density(results$trials.diff.euler), sub=round(p.euler, 5), xlab="Diff.Euler", main= "Euler Characteristic");abline(v=results$real.diff.euler, col="red")


dev.off()

################################### Supplementary Figure 1
pdf("RESULTS/FINAL/SuppFigure_byGene.pdf", height=7.5, width=15)
par(mfrow=c(2,6))

load("RESULTS/FINAL/results_Array_byGene_1000.rda")

results=results_Array_byGene
n_trials=length(results$trials.diff.sdt)
p.sdt=length(which(results$trials.diff.sdt >results$real.diff.sdt))/n_trials
p.euler=length(which(results$trials.diff.euler >results$real.diff.euler))/n_trials

plotPersistenceDiagram(results$ASD.persistance.real, 3,max(results$ASD.persistance.real[,3]), title="Microarray ASD")
plotPersistenceDiagram(results$health.persistance.real, 3, max(results$health.persistance.real[,3]),title="Microrray Controls")
plotPersistenceDiagram(results$ASD.persistance.real, 3,max(results$ASD.persistance.real[which(results$ASD.persistance.real[,1]>0),3]),title="Microrray ASD")
plotPersistenceDiagram(results$health.persistance.real, 3,max(results$health.persistance.real[which(results$health.persistance.real[,1]>0),3]),title="Microrray Controls")
plot(density(results$trials.diff.sdt), sub=round(p.sdt, 5), xlab="Diff.SDT", main= "SDT");abline(v=results$real.diff.sdt, col="red")
plot(density(results$trials.diff.euler), sub=round(p.euler, 5), xlab="Diff.Euler", main= "Euler Characteristic");abline(v=results$real.diff.euler, col="red")

load("RESULTS/FINAL/results_RNAseq_byGene_1000.rda")

results=results_RNAseq_byGene
n_trials=length(results$trials.diff.sdt)
p.sdt=length(which(results$trials.diff.sdt >results$real.diff.sdt))/n_trials
p.euler=length(which(results$trials.diff.euler >results$real.diff.euler))/n_trials

plotPersistenceDiagram(results$ASD.persistance.real, 3,max(results$ASD.persistance.real[,3]), title="RNA-seq ASD")
plotPersistenceDiagram(results$health.persistance.real, 3, max(results$health.persistance.real[,3]),title="RNA-seq Controls")
plotPersistenceDiagram(results$ASD.persistance.real, 3,max(results$ASD.persistance.real[which(results$ASD.persistance.real[,1]>0),3]),title="RNA-seq ASD")
plotPersistenceDiagram(results$health.persistance.real, 3,max(results$health.persistance.real[which(results$health.persistance.real[,1]>0),3]),title="RNA-seq Controls")
plot(density(results$trials.diff.sdt), sub=round(p.sdt, 5), xlab="Diff.SDT", main= "SDT");abline(v=results$real.diff.sdt, col="red")
plot(density(results$trials.diff.euler), sub=round(p.euler, 5), xlab="Diff.Euler", main= "Euler Characteristic");abline(v=results$real.diff.euler, col="red")

dev.off()