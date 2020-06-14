################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### mothur pipeline OTU: rarefaction curves
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2020
### korn@cumulonimbus.at
################################################################################
################################################################################


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("vegan")


rm(list = ls())


setwd("~/Documents/Metabarcoding-Pipeline-Building/")


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
  ### Specifiy path for pro- and eukaryotes
  ifelse(primer == "16S",
         goto <- "prok",
         goto <- "euk")


  ### Read OTU
  if(method == "OTU" && classifier == "RDP")
    {
    wine <- import_mothur(mothur_list_file = NULL,
                          mothur_group_file = NULL,
                          mothur_tree_file = NULL,
                          cutoff = NULL,
                          mothur_shared_file = paste("OTU", goto,
                                                     "wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared",
                                                     sep = "/"),
                          paste("OTU", goto,
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
    seqtab.nochim <- readRDS(paste("ASV", goto, "seq.tab.nochim.rds",
                                   sep = "/"))

    ifelse(classifier == "RDP",
           taxa <- readRDS(paste("ASV", goto, "taxa.rds", sep = "/")),
           taxa <- readRDS(paste("ASV", goto, "taxa.id.rds", sep = "/")))

    wine <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                     tax_table(taxa))

    dna <- Biostrings::DNAStringSet(taxa_names(wine))
    names(dna) <- taxa_names(wine)
    wine <- merge_phyloseq(wine, dna)
    taxa_names(wine) <- paste0("ASV", seq(ntaxa(wine)))

    ### Rename taxonomic ranks
    # colnames(tax_table(wine)) <- c("Domain", "Phylum", "Class", "Order",
    #                               "Family", "Genus")
    colnames(tax_table(wine)) <- c("Domain", "Phylum", "Class", "Order",
                                   "Family", "Genus", "Species")
  }

  ### Transpose (sometimes the OTU table is transposed... :-|)
  if (taxa_are_rows(wine))
    {otu_table(wine) <- t(otu_table(wine))}

  return(wine)
  }


################################################################################
### Read taxa into phyloseq object
## Note: for now choose RDP for OTU, DECIPHER for ASV
method <- "ASV"
primer <- "18S"
classifier <- "RDP"
classifier <- "DECIPHER"


wine <- readTaxa(method = method, primer = primer, classifier = classifier)
wine


colnames(tax_table(wine))


################################################################################
### Add metadata
## To do: merge with metadata
wine <- merge_phyloseq(wine, MetaData)


################################################################################
### Remove spurious taxa
if(primer == "16S")
  {
  wine.s <- subset_taxa(wine, !(Domain %in% c("unknown") |
                                  Phylum %in% c("Eukaryota_unclassified", NA) |
                                  Order %in% c("Chloroplast") |
                                  Family %in% c("Mitochondria")))
  }

if(primer == "18S")
  {
  wine.s <- subset_taxa(wine, !(Domain %in% c("Bacteria", "unknown") |
                                  Phylum %in% c("Eukaryota_unclassified",
                                                "Mollusca", "Vertebrata", NA) |
                                  Class %in% c("Insecta", "Ellipura",
                                               "Embryophyta", "Arachnida",
                                               "Heterophyidae",
                                               # "Ichthyophonae",
                                               "Arthropoda_unclassified",
                                               "unclassified_Hexapoda")))
  }


################################################################################
### Check for empty taxa and remove if any
if(any(taxa_sums(wine) == 0))
  {
  sum(taxa_sums(wine) == 0)
  wine <- prune_taxa(taxa_sums(wine) > 0, wine)
  }


################################################################################
### Percentual abundance
wine.s <- transform_sample_counts(wine.s, function(otu) {otu / sum(otu)})
plot(rowSums(otu_table(wine.s)), ylim = c(0, 1),
     xlab = "Samples", ylab = "Abundance [%]") # 100 % in all samples


################################################################################
### Abundance filtering
wine.a <- filter_taxa(wine.s, function(otu) {mean(otu) > 0.0001},
                      prune = TRUE)
points(rowSums(otu_table(wine.a)), col = "red")


if(any(taxa_sums(wine) == 0))
  {
  sum(taxa_sums(wine) == 0)
  wine <- prune_taxa(taxa_sums(wine) > 0, wine)
  }


################################################################################
### The effect of pruning and filtering
nr.taxa <- data.frame(NrTaxa = c(dim(tax_table(wine))[1],
                                 dim(tax_table(wine.s))[1],
                                 dim(tax_table(wine.a))[1]),
                      Dataset = c("Raw", "Taxonomic filtering",
                                  "Abundance filtering"))
nr.taxa$Dataset <- ordered(nr.taxa$Dataset, levels = c("Raw",
                                                       "Taxonomic filtering",
                                                       "Abundance filtering"))


## To do: change colors
ggplot(nr.taxa, aes(x = Dataset, y = NrTaxa, color = Dataset)) +
  geom_point() +
  theme(legend.position = "none") +
  xlab(expression(paste(alpha, "-Diversity"))) +
  ylab("")


################################################################################
### Plot taxa
plot_bar(wine.a, x = "Phylum", fill = "Phylum") +
  facet_grid(Domain ~ ., scales = "free_y") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity",
           position = "stack") +
  theme(legend.position = "top") +
  xlab("") +
  ylab("Abundance [%]") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  coord_flip()


################################################################################
### Shannon entropy
## To do: add useful metadata
plot_richness(wine.a, x = "Treatment", measures = c("Shannon"),
              color = "Site") +
  geom_boxplot(color = "gray", alpha = 0.1, outlier.shape = NA) +
  facet_wrap( ~ Site, nrow = 1) +
  xlab("") +
  ylab("Shannon entropy") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


################################################################################
### Ordination
## First, log-transform
wine.log <- transform_sample_counts(wine.a, function(otu) {log1p(otu)})

## Nonmetric Multidimensional Scaling (nMDS)
wine.nmds <- ordinate(wine.log, method = "NMDS", distance = "bray", k = 3,
                     autotransform = FALSE, trymax = 100)
wine.nmds

stressplot(wine.nmds)


plot_ordination(wine.a, wine.nmds,# shape = "Succession",
                # color = "Site", title = NULL, label = "Sector", axes = 1:2
                ) +
  # stat_ellipse(aes(group = Site), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal")


## MDS = PCoA
wine.nmds <- ordinate(wine.log, method = "PCoA", distance = "bray")

barplot(wine.nmds$values$Relative_eig)
biplot(wine.nmds, data.frame(otu_table(wine.a)))


plot_ordination(wine.a, wine.nmds,# shape = "Succession",
                # color = "Site", title = NULL, label = "Sector", axes = 1:2
                ) +
  # stat_ellipse(aes(group = Site), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal")


################################################################################
################################################################################
