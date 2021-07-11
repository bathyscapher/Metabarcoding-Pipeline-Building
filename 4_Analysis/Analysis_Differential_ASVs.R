###############################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### Analysis of microbiome data
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2021
### korn@cumulonimbus.at
################################################################################
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20))


## Install (if necessary) and load the DeSeq2 library
suppressMessages(deseq2_installed <- require(DESeq2))
if (!deseq2_installed) {
  BiocManager::install("DESeq2")
  }

library("DESeq2")
packageVersion("DESeq2")


rm(list = ls())
setwd("/home/rstudio/")


################################################################################
### The readTaxa function expects the data to be arranged in the directory as:
### wd: the given working directory
### |-- OTU
###    |-- prok: 16S samples
###    |-- euk: 18S samples
### |-- ASV
###    |-- prok: 16S samples
###    |-- euk: 18S samples


readTaxa <- function(method = c("OTU", "ASV"), primer = c("16S", "18S"),
                     classifier = c("RDP", "DECIPHER")){
  ### Specify path for pro- and eukaryotes
  ifelse(primer == "16S",
         goto <- "prok/filtered",
         goto <- "euk/filtered")
  
  
  ### Read OTU
  if(method == "OTU" && classifier == "RDP")
  {
    wine <- import_mothur(mothur_list_file = NULL,
                          mothur_group_file = NULL,
                          mothur_tree_file = NULL,
                          cutoff = NULL,
                          mothur_shared_file = paste(goto,
                                                     "wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared",
                                                     sep = "/"),
                          paste(goto,
                                "wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy",
                                sep = "/"),
                          parseFunction = parse_taxonomy_default)
    
    ## Rename taxonomic ranks
    colnames(tax_table(wine)) <- c("Domain", "Phylum", "Class", "Order",
                                   "Family", "Genus")
  }
  
  
  ### Read ASV
  if (method == "ASV")
  {
    seqtab.nochim <- readRDS(paste(goto, "asv.tab.nochim.rds",
                                   sep = "/"))
    
    ifelse(classifier == "RDP",
           taxa <- readRDS(paste(goto, "taxa.rds", sep = "/")),
           taxa <- readRDS(paste(goto, "taxa.id.rds", sep = "/")))
    
    wine <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                     tax_table(taxa))
    
    dna <- Biostrings::DNAStringSet(taxa_names(wine))
    names(dna) <- taxa_names(wine)
    wine <- merge_phyloseq(wine, dna)
    taxa_names(wine) <- paste0("ASV", seq(ntaxa(wine)))
    
    ### Rename taxonomic ranks
    ifelse(classifier == "RDP",
           colnames(tax_table(wine)) <- c("Domain", "Phylum", "Class", "Order",
                                          "Family", "Genus"),
           colnames(tax_table(wine)) <- c("Domain", "Phylum", "Class", "Order",
                                          "Family", "Genus", "Species"))
  }
  
  ### Transpose (sometimes the OTU table is transposed... :-|)
  if (taxa_are_rows(wine))
  {otu_table(wine) <- t(otu_table(wine))}
  
  return(wine)
}


## Read the metadata of the study 
metadata <- read.csv("Metabarcoding-Pipeline-Building/Metadata_corr.csv",
                     sep = ",", header = TRUE, row.names = 1)
head(metadata)


metadata$treatment[metadata$treatment == 'AC'] <- 'AlternatingCover'
metadata$treatment[metadata$treatment == 'BG'] <- 'BareGround'
metadata$treatment[metadata$treatment == 'CC'] <- 'CompleteCover'


## Code vineyard names as factors
metadata$treatment <- as.factor(metadata$treatment)


## Rename and sort factors
metadata$ord.treatment <- factor(metadata$treatment)
levels(metadata$ord.treatment) <- c("AlternatingCover", "BareGround",
                                    "CompleteCover")
metadata$ord.treatment <- ordered(metadata$ord.treatment,
                                  levels = c("BareGround", "AlternatingCover",
                                             "CompleteCover"))


## Read OTU/ASV table
tax.asv <- readTaxa("ASV", "16S", "RDP")
sample_data(tax.asv) <- metadata
tax.asv


tax.otu <- readTaxa("OTU", "16S", "RDP")
tax.otu <- prune_taxa(taxa_sums(tax.otu) > 1, tax.otu)
sample_data(tax.otu) <- metadata
tax.otu


## Determine differentially abundant taxa based on the treatment
## Note: DESeq2 does not handle ordered factors

tax.asv.deseq <- phyloseq_to_deseq2(tax.asv, ~treatment)


## Function to calculate the geometric mean
gm_mean <- function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }


## Calculate the geometric mean
geoMeans <- apply(counts(tax.asv.deseq), 1, gm_mean)
tax.asv.deseq.sf <- estimateSizeFactors(tax.asv.deseq, geoMeans = geoMeans)

tax.asv.deseq.de <- DESeq(tax.asv.deseq.sf, fitType = "local")


## Order results by decreasing adjusted p-value
res <- results(tax.asv.deseq.de)
res <- res[order(res$padj, na.last = NA), ]
dim(res)
head(res)


## Extract all taxa differentially abundant at a significance level 0.05
alpha <- 0.05
sigtab.raw <- res[(res$padj < alpha), ]
dim(sigtab.raw)
sigtab <- cbind(as(sigtab.raw, "data.frame"),
                as(tax_table(tax.asv)[rownames(sigtab.raw), ], "matrix"))
head(sigtab)
dim(sigtab)


## Plot the fold-change of the DA taxa with the genus information and color
## by phylum
sigtabgen <- subset(sigtab, !is.na(Genus))
dim(sigtabgen)


## Phylum order: sort by log2FoldChange and sort phyla accordingly
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x <- sort(x, TRUE)


sigtabgen$Phylum <- factor(as.character(sigtabgen$Phylum), levels = names(x))


## Genus order: sort by log2FoldChange and sort genera accordingly
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x <- sort(x, TRUE)


sigtabgen$Genus <- factor(as.character(sigtabgen$Genus), levels = names(x))


## Plot
ggplot(sigtabgen, aes(y = Genus, x = log2FoldChange, color = Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size = 6) +
  labs(title = "CompleteCover vs AlternatingCover") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))


################################################################################
## Exercise: DA between 2 groups only
table(sample_data(tax.asv)$treatment)
subset_samples(tax.asv, sample_data(tax.asv)$treatment != "CompleteCover")


################################################################################
################################################################################