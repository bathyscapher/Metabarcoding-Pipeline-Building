################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### dada2 pipeline
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2021
### korn@cumulonimbus.at
################################################################################
################################################################################


library("dada2")
library("DECIPHER")


rm(list = ls())


## Choose pro- or eukaryotes
# primer <- "16S"
primer <- "18S"


if (primer == "16S") {
  setwd("/home/rstudio/prok/filtered/")
  } else {
    setwd("/home/rstudio/euk/filtered")
    }

getwd()


################################################################################
ncore <- 6 # number of available processors

list.files(pattern = "fastq.gz")


rF <- sort(list.files(pattern = "_R1_filt.fastq.gz", full.names = TRUE))
rR <- sort(list.files(pattern = "_R2_filt.fastq.gz", full.names = TRUE))


sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 1)


names(rF) <- sample.names
names(rR) <- sample.names


### Estimate and plot the error rates
errF <- learnErrors(rF, multithread = ncore, verbose = TRUE)
errR <- learnErrors(rR, multithread = ncore, verbose = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)


### Core sample inference algorithm
dadaF <- dada(rF, err = errF, multithread = ncore, verbose = TRUE)
dadaR <- dada(rR, err = errR, multithread = ncore, verbose = TRUE)


saveRDS(dadaF, "dadaF.rds")
saveRDS(dadaR, "dadaR.rds")


### Merge paired reads
contigs <- mergePairs(dadaF, rF, dadaR, rR, verbose = TRUE)
head(contigs[[1]])

saveRDS(contigs, "contigs.rds")


### Construct ASV table
asv.tab <- makeSequenceTable(contigs)
dim(asv.tab)

table(nchar(getSequences(asv.tab)))


### Chimera detection
asv.tab.nochim <- removeBimeraDenovo(asv.tab, method = "consensus",
                                     multithread = ncore, verbose = TRUE)
dim(asv.tab.nochim)
sum(asv.tab.nochim) / sum(asv.tab)


table(nchar(getSequences(asv.tab.nochim)))


saveRDS(asv.tab.nochim, "asv.tab.nochim.rds")


### Assign taxonomy with RDP classifier
taxa <- assignTaxonomy(asv.tab.nochim,
                       "/home/rstudio/silva_nr99_v138.1_train_set.fa.gz",
                       multithread = TRUE, verbose = TRUE)


saveRDS(taxa, "taxa.rds")


### Assign taxonomy with IdTaxa
dna <- DNAStringSet(getSequences(asv.tab.nochim))

load("/home/rstudio/SILVA_SSU_r138_2019.RData")


ids <- IdTaxa(dna, trainingSet, strand = "top", processors = ncore,
              verbose = TRUE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")


## Convert output object of class "Taxa"
taxa.id <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  # taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))


colnames(taxa.id) <- ranks
rownames(taxa.id) <- getSequences(asv.tab.nochim)


saveRDS(taxa.id, "taxa.id.rds")


################################################################################
################################################################################
