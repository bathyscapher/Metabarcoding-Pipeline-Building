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
library("reshape2")
library("ggplot2")


rm(list = ls())


setwd("/scratch/PromESSinG/prok/filtered/")
# setwd("/scratch/PromESSinG/euk/filtered/")


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


### Core sample inference algorithm
dadaF <- dada(rF, err = errF, multithread = ncore)
dadaR <- dada(rR, err = errR, multithread = ncore)


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
                       "~/Desktop/SMP_unsynced/silva/silva_nr_v138_train_set.fa.gz",
                       multithread = TRUE, verbose = TRUE)


saveRDS(taxa, "taxa.rds")


### Assign taxonomy with IdTaxa
dna <- DNAStringSet(getSequences(asv.tab.nochim))

load("~/Desktop/SMP_unsynced/silva/SILVA_SSU_r132_March2018.RData")


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


### Survey where the reads are 'lost' in the pipeline.
getN <- function(x) {
  sum(getUniques(x))
}


setwd("/scratch/PromESSinG/")


## 16S
FASTQ.f <- readRDS("prok/FASTQ.f.rds")
dadaF <- readRDS("prok/filtered/dadaF.rds")
dadaR <- readRDS("prok/filtered/dadaR.rds")
contigs <- readRDS("prok/filtered/contigs.rds")
asv.tab.nochim <- readRDS("prok/filtered/asv.tab.nochim.rds")
# taxa <- readRDS("prok/filtered/taxa.rds")
# taxa.id <- readRDS("prok/filtered/taxa.id.rds")

track.prok <- cbind(FASTQ.f, sapply(dadaF, getN), sapply(dadaR, getN),
                   sapply(contigs, getN), rowSums(asv.tab.nochim))

colnames(track.prok) <- c("Raw", "Filtered", "Denoised F", "Denoised R",
                         "Contigs", "Without Chimeras")
rownames(track.prok) <- sample.names

track.cS.prok <- as.data.frame(t(colSums(track.prok)))
track.cS.prok$Primer <- "16S"


## 18S
FASTQ.f <- readRDS("euk/FASTQ.f.rds")
dadaF <- readRDS("euk/filtered/dadaF.rds")
dadaR <- readRDS("euk/filtered/dadaR.rds")
contigs <- readRDS("euk/filtered/contigs.rds")
asv.tab.nochim <- readRDS("euk/filtered/asv.tab.nochim.rds")
# taxa <- readRDS("euk/filtered/taxa.rds")
# taxa.id <- readRDS("euk/filtered/taxa.id.rds")


track.euk <- cbind(FASTQ.f, sapply(dadaF, getN), sapply(dadaR, getN),
                   sapply(contigs, getN), rowSums(asv.tab.nochim))

colnames(track.euk) <- c("Raw", "Filtered", "Denoised F", "Denoised R",
                         "Contigs", "Without Chimeras")
rownames(track.euk) <- sample.names

track.cS.euk <- as.data.frame(t(colSums(track.euk)))
track.cS.euk$Primer <- "18S"


## Combine
track.cS <- rbind(track.cS.prok, track.cS.euk)
track.cS.m <- melt(track.cS, variable.name = "dada2", value.name = "Reads")


## Plot
ggplot(data = track.cS.m, aes(x = dada2, y = log(Reads), color = Primer)) +
  geom_line(aes(group = Primer), linetype = "dashed") +
  geom_point() +
  facet_grid( ~ Primer) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("")


################################################################################
################################################################################
