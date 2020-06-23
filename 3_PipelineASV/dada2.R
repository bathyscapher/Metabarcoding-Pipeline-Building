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


### Core sample inference algorithm
dadaF <- dada(rF.f, err = errF, multithread = ncore)
dadaR <- dada(rR.f, err = errR, multithread = ncore)


saveRDS(dadaF, "dadaF.rds")
saveRDS(dadaR, "dadaR.rds")


### Merge paired reads
contigs <- mergePairs(dadaF, rF.f, dadaR, rR.f, verbose = TRUE)
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

saveRDS(asv.tab.nochim, "asv.tab.nochim.rds")


### Assign taxonomy with IdTaxa and SILVA
dna <- DNAStringSet(getSequences(asv.tab.nochim))

load("SILVA_SSU_r132_March2018.RData")


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


setwd("...")

## 16S
FASTQ.f <- readRDS("16S/FASTQ.f.rds")
dadaF <- readRDS("16S/filtered/dadaF.rds")
dadaR <- readRDS("16S/filtered/dadaR.rds")
contigs <- readRDS("16S/filtered/contigs.rds")
asv.tab.nochim <- readRDS("16S/filtered/seq.tab.nochim.rds")
taxa.id <- readRDS("16S/filtered/taxa.id.rds")

track.16S <- cbind(FASTQ.f, sapply(dadaF, getN), sapply(dadaR, getN),
                   sapply(contigs, getN), rowSums(asv.tab.nochim))

colnames(track.16S) <- c("Raw", "Filtered", "Denoised F", "Denoised R",
                         "Contigs", "Without Chimeras")
rownames(track.16S) <- sample.names

track.cS.16S <- as.data.frame(t(colSums(track.16S)))
track.cS.16S$Primer <- "16S"


## 18S
FASTQ.f <- readRDS("18S/FASTQ.f.rds")
dadaF <- readRDS("18S/filtered/dadaF.rds")
dadaR <- readRDS("18S/filtered/dadaR.rds")
contigs <- readRDS("18S/filtered/contigs.rds")
asv.tab.nochim <- readRDS("18S/filtered/seq.tab.nochim.rds")
taxa.id <- readRDS("18S/filtered/taxa.id.rds")

track.18S <- cbind(FASTQ.f, sapply(dadaF, getN), sapply(dadaR, getN),
                   sapply(contigs, getN), rowSums(asv.tab.nochim))

colnames(track.18S) <- c("Raw", "Filtered", "Denoised F", "Denoised R",
                         "Contigs", "Without Chimeras")
rownames(track.18S) <- sample.names

track.cS.18S <- as.data.frame(t(colSums(track.18S)))
track.cS.18S$Primer <- "18S"


## Combine
track.cS <- rbind(track.cS.16S, track.cS.18S)
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
