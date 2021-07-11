###############################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### OTUs vs. ASVs
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2021
### korn@cumulonimbus.at
################################################################################
################################################################################


library("phyloseq")
library("VennDiagram")
library("gridExtra")


rm(list = ls())
setwd("/home/rstudio/")


## Colorblind scale
cols <- c("#C59434", "#999999", "#009E73") # brown, gray, green


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


################################################################################
### Read taxa into phyloseq object ####
tax.asv <- readTaxa("ASV","16S","RDP")
tax.otu <- readTaxa("OTU","16S","RDP")


################################################################################
### Remove singletons ###
sum(taxa_sums(tax.asv) == 1)
# tax.asv <- prune_taxa(taxa_sums(tax.asv) > 1, tax.asv)
# tax.asv


sum(taxa_sums(tax.otu) == 1)
tax.otu <- prune_taxa(taxa_sums(tax.otu) > 1, tax.otu)
tax.otu


################################################################################
tax.asv.tax <- as.data.frame(tax_table(tax.asv))
tax.otu.tax <- as.data.frame(tax_table(tax.otu))


## On the genus level
s1 <- unique(tax.asv.tax$Genus)
s2 <- unique(tax.otu.tax$Genus)
over <- intersect(s1,s2)


dev.off()
venn <- draw.pairwise.venn(length(s1), length(s2), length(over), fill=c(2,3),
                           cat.col=c(2,3), category=c("ASV", "OTU"),)
grid.arrange(gTree(children = venn), top = "On the genus-level")


## On the family level
s1 <- unique(tax.asv.tax$Family)
s2 <- unique(tax.otu.tax$Family)
over <- intersect(s1,s2)


dev.off()
venn <- draw.pairwise.venn(length(s1), length(s2), length(over), fill=c(2,3),
                           cat.col=c(2,3), category=c("ASV", "OTU"))
grid.arrange(gTree(children = venn), top = "On the family-level")


## On the order level
s1 <- unique(tax.asv.tax$Order)
s2 <- unique(tax.otu.tax$Order)

over <- intersect(s1,s2)

dev.off()
venn <- draw.pairwise.venn(length(s1), length(s2), length(over), fill=c(2,3),
                           cat.col=c(2,3), category=c("ASV", "OTU"))
grid.arrange(gTree(children = venn), top = "On the order-level")


################################################################################
### Add metadata
metadata <- read.csv("Metabarcoding-Pipeline-Building/Metadata_corr.csv",
                     sep = ",", header = TRUE, row.names = 1)

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


tax.otu <- merge_phyloseq(tax.otu, sample_data(metadata))
tax.asv <- merge_phyloseq(tax.asv, sample_data(metadata))


### Compare OTU and ASV communities in a nMDS
set.seed(11948)

otu.nmds <- ordinate(tax.otu, "NMDS", "bray", autotransform = FALSE,
                     trymax = 50, k = 2)
otu.nmds
stressplot(otu.nmds)


set.seed(11948)
asv.nmds <- ordinate(tax.asv, "NMDS", "bray", autotransform = FALSE,
                     trymax = 50, k = 2)
asv.nmds
stressplot(asv.nmds)


p.asv.nmds <- plot_ordination(tax.asv, asv.nmds, color = "ord.treatment",
                shape = "ord.treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  scale_color_manual(values = cols) +
  geom_text(aes(label = sample_data(tax.asv)$vineyard), color = "black",
            size = 3) +
  ggtitle("NMDS of Bray-Curtis distance") +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = ord.treatment), type = "t", linetype = 2, size = 0.2)


p.otu.nmds <- plot_ordination(tax.otu, otu.nmds, color = "ord.treatment",
                shape = "ord.treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  scale_color_manual(values = cols) +
  geom_text(aes(label = sample_data(tax.otu)$vineyard), color = "black",
            size = 3) +
  ggtitle("NMDS of Bray-Curtis distance") +
  theme(legend.position = "none") +
  
  coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = ord.treatment), type = "t", linetype = 2, size = 0.2)


grid.arrange(p.otu.nmds, p.asv.nmds, ncol = 2)


################################################################################
################################################################################
