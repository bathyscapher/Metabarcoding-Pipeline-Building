# dada2 pipeline
Clear workspace, set working directory and specify number of available processors.
```R
rm(list = ls())
setwd("...")

library("dada2")
library("DECIPHER")

ncore <- 6 # number of available cores
```

## List fastq files
List all *.fastq.gz files in the working directory and get the sample names.
```R
list.files(pattern = "fastq.gz")

rF <- sort(list.files(pattern = "_R1.fastq.gz", full.names = TRUE))
rR <- sort(list.files(pattern = "_R2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(rF), "_"), `[`, 1)

rF.f <- file.path("filtered", paste0(sample.names, "_16S_R1_filt.fastq.gz"))
rR.f <- file.path("filtered", paste0(sample.names, "_16S_R2_filt.fastq.gz"))

names(rF.f) <- sample.names
names(rR.f) <- sample.names
```

## Error rates
Estimate and plot the error rates.
```
errF <- learnErrors(rF.f, multithread = ncore)
errR <- learnErrors(rR.f, multithread = ncore)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
```

![dada2 error rates 16S](/Graphs/dada2_ErrorRates_16S.png)
![dada2 error rates 18S](/Graphs/dada2_ErrorRates_18S.png)

## Core sample inference algorithm
```R
dadaF <- dada(rF.f, err = errF, multithread = ncore)
dadaR <- dada(rR.f, err = errR, multithread = ncore)
```

Save (and read) intermediate results.
```R
saveRDS(dadaF, "dadaF.rds")
saveRDS(dadaR, "dadaR.rds")
```

## Merge paired reads
```R
contigs <- mergePairs(dadaF, rF.f, dadaR, rR.f, verbose = TRUE)
head(contigs[[1]])

saveRDS(contigs, "contigs.rds")
```

## Construct ASV table
```R
asv.tab <- makeSequenceTable(contigs)
dim(asv.tab)
```

Show distribution of sequence lengths.
```R
table(nchar(getSequences(asv.tab)))
```

## Chimera detection
```R
asv.tab.nochim <- removeBimeraDenovo(asv.tab, method = "consensus",
                                     multithread = ncore, verbose = TRUE)
dim(asv.tab.nochim)

sum(asv.tab.nochim) / sum(asv.tab)

saveRDS(asv.tab.nochim, "asv.tab.nochim.rds")
```

## Assign taxonomy with IdTaxa and SILVA
Convert the chimera-cleaned sequences into a 'DNAStringSet', load the SILVA db and classify the sequences.
```R
dna <- DNAStringSet(getSequences(asv.tab.nochim))

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
rownames(taxid) <- getSequences(asv.tab.nochim)

saveRDS(taxid, "taxaid.rds")
```

## Track reads
Survey where the reads are 'lost' in the pipeline.
```R
getN <- function(x){
  sum(getUniques(x))
  }

track <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN),
               sapply(contigs, getN), rowSums(asv.tab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)
```

![dada2 track reads](/Graphs/dada2_TrackReads.png)
