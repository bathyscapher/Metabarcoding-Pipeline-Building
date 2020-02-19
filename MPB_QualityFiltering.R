################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### Quality filtering with dada2
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2020
### korn@cumulonimbus.at
################################################################################
################################################################################

rm(list = ls())

setwd("/scratch/PromESSinG/16S/")


library("dada2") #; packageVersion("dada2")
library("DECIPHER")#; packageVersion("DECIPHER")
library("gridExtra")


################################################################################
ncore <- 9 # number of available processors


list.files(pattern = "fastq.gz")

rF <- sort(list.files(pattern = "_R1_sub.fastq.gz", full.names = TRUE))
rR <- sort(list.files(pattern = "_R2_sub.fastq.gz", full.names = TRUE))

## Extract sample names
sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 1)

## Filter and trim. Place filtered files in "filtered" subdirectory.
rF.f <- file.path("filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
rR.f <- file.path("filtered", paste0(sample.names, "_R2_filt.fastq.gz"))

names(rF.f) <- sample.names
names(rR.f) <- sample.names


## Filter
FASTQ.f <- filterAndTrim(rF, rF.f, rR, rR.f,
                         maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                         compress = TRUE, multithread = ncore, verbose = TRUE)


### Plot quality profiles
set.seed(74398)
startPlot <- ceiling(runif(1, 0, length(rF) - 3))


plot.rF <- plotQualityProfile(rF[startPlot:(startPlot + 2)])
plot.rF.f <- plotQualityProfile(rF.f[startPlot:(startPlot + 2)])

plot.rR <- plotQualityProfile(rR[startPlot:(startPlot + 2)])
plot.rR.f <- plotQualityProfile(rR.f[startPlot:(startPlot + 2)])

grid.arrange(plot.rF, plot.rR, plot.rF.f, plot.rR.f, ncol = 2)


################################################################################
################################################################################
