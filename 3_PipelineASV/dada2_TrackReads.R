################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### Quality filtering with dada2
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2020
### korn@cumulonimbus.at
################################################################################
################################################################################


library("reshape2")
library("ggplot2")


rm(list = ls())


setwd("/scratch/mpb/")


################################################################################
### Survey where the reads are 'lost' in the pipeline.
getN <- function(x) {
  sum(getUniques(x))
}


rF <- sort(list.files(path = "prok", pattern = "_R1_filt.fastq.gz",
                      full.names = TRUE))
sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 1)


## 16S
FASTQ.f <- readRDS("prok/FASTQ.f.rds")
dadaF <- readRDS("prok/filtered/dadaF.rds")
dadaR <- readRDS("prok/filtered/dadaR.rds")
contigs <- readRDS("prok/filtered/contigs.rds")
asv.tab.nochim <- readRDS("prok/filtered/asv.tab.nochim.rds")


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
        # legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  theme_bw()


################################################################################
################################################################################
