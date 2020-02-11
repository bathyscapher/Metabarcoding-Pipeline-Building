# dada2 pipeline
Clear workspace and set working directory and specify number of available processors.
```R
rm(list = ls())
setwd("...")

library("dada2")
library("DECIPHER")
library("ggplot2")
  theme_set(theme_bw(base_size = 15) +
    theme(rect = element_rect(fill = "transparent")))
library("phyloseq")

ncore <- 6 # specify number of available cores
```

## Read fastq files
List all *.fastq.gz files in the working directory and get the sample names.
```R
list.files(pattern = "fastq.gz")

rF <- sort(list.files(pattern = "_R1.fastq.gz", full.names = TRUE))
rR <- sort(list.files(pattern = "_R2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 1)
```


## Quality control
Filter  and trim the reads. Place filtered files in filtered/ subdirectory. Compare read quality profiles before and after filtering exemplarily.
```R
rF.f <- file.path("filtered", paste0(sample.names, "_16S_R1_filt.fastq.gz"))
rR.f <- file.path("filtered", paste0(sample.names, "_16S_R2_filt.fastq.gz"))

names(rF.f) <- sample.names
names(rR.f) <- sample.names

out <- filterAndTrim(rF, rF.f, rR, rR.f,
                     maxN = 0, maxEE = c(2,2), truncQ = 10, rm.phix = TRUE,
                     compress = TRUE, multithread = ncore)

plotQualityProfile(rF[1:3])
plotQualityProfile(rF.f[1:3])
plotQualityProfile(rR[1:3])
plotQualityProfile(rR.f[1:3])
```


## Error rates
Estimate and plot the error rates.
```
errF <- learnErrors(rF.f, multithread = ncore)
errR <- learnErrors(rR.f, multithread = ncore)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
```


## Core sample inference algorithm
```R
dadaF <- dada(rF.f, err = errF, multithread = ncore)
dadaR <- dada(rR.f, err = errR, multithread = ncore)
# dadaFs[[1]]
```

Save (and read) intermediate results.
```R
saveRDS(dadaF, "dadaF.rds")
saveRDS(dadaR, "dadaR.rds")
# dadaF <- readRDS("dadaF.rds")
# dadaR <- readRDS("dadaR.rds")
```


## Merge paired reads
```R
contigs <- mergePairs(dadaF, rF.f, dadaR, rR.f, verbose = TRUE)
head(contigs[[1]])

saveRDS(contigs, "contigs.rds")
# contigs <- readRDS("contigs.rds")
```

## Construct ASV table
```R
seqtab <- makeSequenceTable(contigs)
dim(seqtab)
```

Show distribution of sequence lengths.
```R
table(nchar(getSequences(seqtab)))
```

## Chimera detection
```R
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = ncore, verbose = TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim) / sum(seqtab)

saveRDS(seqtab.nochim, "seqtab.nochim.rds")
# seqtab.nochim <- readRDS("seqtab.nochim.rds")
```

## Track reads through the pipeline
To survey where the reads are 'lost'.
```R
getN <- function(x){
  sum(getUniques(x))
  }

track <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN),
               sapply(contigs, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)
``````

## Assign taxonomy with IdTaxa and SILVA
Convert the chimera-cleaned sequences into a 'DNAStringSet', load the SILVA db and classify the sequences.
```R
dna <- DNAStringSet(getSequences(seqtab.nochim))

load("SILVA_SSU_r132_March2018.RData")

ids <- IdTaxa(dna, trainingSet, strand = "top", processors = ncore,
              verbose = TRUE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

## Convert output object of class "Taxa"
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  # taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
  }))

colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)

saveRDS(taxid, "taxaid.rds")
# taxaid <- readRDS("taxaid.rds")
```
