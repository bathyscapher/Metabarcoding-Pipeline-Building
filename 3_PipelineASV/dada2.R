################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### dada2 pipeline
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2021
### korn@cumulonimbus.at
################################################################################
################################################################################


library("dada2")
packageVersion("dada2")
library("DECIPHER")


rm(list = ls())


## Choose pro- or eukaryotes
primer <- "16S"
# primer <- "18S"


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
sample.names


names(rF) <- sample.names
names(rR) <- sample.names


### Estimate and plot the error rates
## Increase number of bases used to learn the errors increases model accuracy
errF <- learnErrors(rF, multithread = ncore, verbose = TRUE, nbases = 1e10)
errR <- learnErrors(rR, multithread = ncore, verbose = TRUE, nbases = 1e10)


saveRDS(errF, "errF.rds")
saveRDS(errR, "errR.rds")
# errF <- readRDS("errF.rds")
# errR <- readRDS("errR.rds")


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
asv.tab.nochim <- readRDS("asv.tab.nochim.rds")


taxa <- assignTaxonomy(asv.tab.nochim, 
                       # "/mothur/refs/silva_nr_v132_train_set.fa.gz", # 18S
                       "/mothur/refs/silva_nr99_v138.1_train_set.fa.gz", # 16S
                       multithread = TRUE, verbose = TRUE)


saveRDS(taxa, "taxa.rds")


### Add species
taxa.species <- assignSpecies(taxa,
                              "/mothur/refs/silva_species_assignment_v138.1.fa.gz",
                              allowMultiple = TRUE, verbose = TRUE)

saveRDS(taxa.species, "taxa.species.rds")


taxa.species <- readRDS("taxa.species.rds")
dim(taxa.species)


### Assign taxonomy with IdTaxa
dna <- DNAStringSet(getSequences(asv.tab.nochim))


load("/mothur/refs/SILVA_SSU_r138_2019.RData")


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
## Exercise: export the first 10 sequences and identify them with BLAST
## https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
names(dna) <- sprintf("ASV%04d", 1:length(dna))
writeXStringSet(dna[1:10], "asv.tab.nochim.10.FASTA")


################################################################################
################################################################################
