################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### Quality filtering with dada2
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2020
### korn@cumulonimbus.at
################################################################################
################################################################################


library("dada2")
library("DECIPHER")


rm(list = ls())


setwd("/scratch/PromESSinG/euk/filtered/")


################################################################################
ncore <- 9 # number of available processors

list.files(pattern = "fastq.gz")


rF <- sort(list.files(pattern = "_R1_filt.fastq.gz", full.names = TRUE))
rR <- sort(list.files(pattern = "_R2_filt.fastq.gz", full.names = TRUE))


sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 1)


names(rF) <- sample.names
names(rR) <- sample.names


### Estimate and plot the error rates

errF <- learnErrors(rF, multithread = ncore)
errR <- learnErrors(rR, multithread = ncore)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)



################################################################################
################################################################################
