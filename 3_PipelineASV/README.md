# Pipeline building ASV
Generally following the [dada2 tutorial](https://benjjneb.github.io/dada2/tutorial.html), the major steps are:

1. Estimate error rates
1. Core sample inference algorithm
1. Merge reads
1. Construct ASV table
1. Remove chimeras
1. Assign taxonomy (RDP classifier and IdTaxa)
1. Track reads


## dada2 pipeline
Clear workspace, set working directory and specify number of available processors.
```R
library("dada2")
library("DECIPHER")

rm(list = ls())

ncore <- 6
```


## Choose pro- or eukaryotes
```
primer <- "16S"
# primer <- "18S"


if (primer == "16S") {
  setwd("/home/rstudio/prok/filtered/")
  } else {
    setwd("/home/rstudio/euk/filtered")
  }
getwd()
```

## List fastq files
List all `*.fastq.gz` files in the working directory and get the sample names.
```R
list.files(pattern = "fastq.gz")

rF <- sort(list.files(pattern = "_R1_filt.fastq.gz", full.names = TRUE))
rR <- sort(list.files(pattern = "_R2_filt.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 1)
sample.names

names(rF.f) <- sample.names
names(rR.f) <- sample.names
```


## Error rates
Estimate and plot the error rates.
```
errF <- learnErrors(rF, multithread = ncore)
errR <- learnErrors(rR, multithread = ncore)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
```

![dada2 error rates 16S](/Graphs/dada2_ErrorRates_16S.png)
![dada2 error rates 18S](/Graphs/dada2_ErrorRates_18S.png)


## Core sample inference algorithm
```R
dadaF <- dada(rF, err = errF, multithread = ncore, verbose = TRUE)
dadaR <- dada(rR, err = errR, multithread = ncore, verbose = TRUE)
```

Save intermediate results.
```R
saveRDS(dadaF, "dadaF.rds")
saveRDS(dadaR, "dadaR.rds")
```


## Merge paired reads
```R
contigs <- mergePairs(dadaF, rF, dadaR, rR, verbose = TRUE)
head(contigs[[1]])

saveRDS(contigs, "contigs.rds")
```


## Construct ASV table
Build ASV table and show dimension and distribution of sequence lengths.
```R
asv.tab <- makeSequenceTable(contigs)
dim(asv.tab)

table(nchar(getSequences(asv.tab)))
```

## Chimera detection
```R
asv.tab.nochim <- removeBimeraDenovo(asv.tab, method = "consensus",
                                     multithread = ncore, verbose = TRUE)
dim(asv.tab.nochim)
sum(asv.tab.nochim) / sum(asv.tab)

table(nchar(getSequences(asv.tab.nochim)))

saveRDS(asv.tab.nochim, "asv.tab.nochim.rds")
```

## Assign taxonomy
Convert the chimera-free sequences into a `DNAStringSet`, load the SILVA db and classify the sequences. There are two classifiers implemented in `dada2`: the RDP classifier and IdTaxa.

### RDP classifier
```R
asv.tab.nochim <- readRDS(file = "asv.tab.nochim.rds")

taxa <- assignTaxonomy(asv.tab.nochim, "silva_nr99_v138.1_train_set.fa.gz",
                       multithread = TRUE, verbose = TRUE)

saveRDS(taxa, "taxa.rds")
```


### Add species
```R
taxa.species <- assignSpecies(taxa,
                              "silva_species_assignment_v138.1.fa.gz",
                              verbose = TRUE)

saveRDS(taxa.species, "taxa.species.rds")
```


### IdTaxa

```R
dna <- DNAStringSet(getSequences(asv.tab.nochim))

load("/home/rstudio/silva/old_SILVA_SSU_r138_2019.RData")


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
```


## Track reads
Survey where the reads are 'lost' in the pipeline.
```R
getN <- function(x) {
  sum(getUniques(x))
  }


setwd("...")

### 16S
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

### 18S
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

### Combine
track.cS <- rbind(track.cS.16S, track.cS.18S)
track.cS.m <- melt(track.cS, variable.name = "dada2", value.name = "Reads")

### Plot
ggplot(data = track.cS.m, aes(x = dada2, y = log(Reads), color = Primer)) +
  geom_line(aes(group = Primer), linetype = "dashed") +
  geom_point() +
  facet_grid( ~ Primer) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("")
```

![dada2 track reads](/Graphs/dada2_TrackReads.png)
