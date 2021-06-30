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
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("vegan")


rm(list = ls())
setwd("~") # set working directory to home = /home/rstudio
setwd("~/docker/")


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
## Note: for now choose RDP for OTU, DECIPHER for ASV
method <- "ASV"
primer <- "18S"
classifier <- "RDP"
# classifier <- "DECIPHER"


wine <- readTaxa(method = method, primer = primer, classifier = classifier)
wine


colnames(tax_table(wine))


################################################################################
### Add metadata ####
## Load metadata: adjust path if necessary, file is online on
## https://github.com/bathyscapher/Metabarcoding-Pipeline-Building
metadata <- read.csv("Metabarcoding-Pipeline-Building/Metadata_corr.csv",
                     sep = ",", header = TRUE, row.names = 1)
head(metadata)


## Code vineyard names as factors
metadata$vineyard <- as.factor(metadata$vineyard)


## Rename and sort factors
metadata$treatment <- factor(metadata$treatment)
levels(metadata$treatment) <- c("Alternating cover", "Bare ground",
                                "Complete cover")
metadata$treatment <- ordered(metadata$treatment,
                              levels = c("Bare ground", "Alternating cover",
                                         "Complete cover"))


MetaData <- sample_data(metadata)
wine <- merge_phyloseq(wine, MetaData)


cols <- c("#C59434", "#999999", "#009E73") # brown, gray, green


### Plot taxa ####
plot_bar(wine, x = "Phylum") +
  facet_grid(Domain ~ ., scales = "free_y", space = "free") +
  geom_bar( stat = "identity",
           position = "stack") +
  theme(legend.position = "top", legend.direction = "horizontal") +
  xlab("") +
  ylab("Abundance [%]") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  coord_flip()



## Exploring our data
nsamples(wine) # number of samples
ntaxa(wine) # number of taxa
sample_names(wine)[1:66]
sample_variables(wine) # metadata variables
sample_data(wine)
otu_table(wine)[1:4, 1:5]
tax_table(wine)[1:5, 1:6]


## Which taxa do we have in our data? Here on the phlyum level:
get_taxa_unique(wine, taxonomic.rank = rank_names(wine)[2], errorIfNULL = TRUE)


## Alternatively, the same and sorted
sort(unique(tax_table(wine)[, 2]))


## Look at our data
readsumsdf <- data.frame(nreads = sort(taxa_sums(wine), TRUE),
                         sorted = 1:ntaxa(wine), type = "OTUs")


## Plot our raw reads per OTU
ggplot(readsumsdf, aes(x = sorted, y = nreads)) +
  geom_bar(stat = "identity") +
  ggtitle("Total number of reads") +
  scale_y_log10() +
  xlab("") +
  ylab(expression(paste(log[10], "(Number of reads)")))


wine
colSums(readsumsdf[, 1, drop = FALSE])


################################################################################
### Taxonomic filtering: remove spurious taxa ####
if(primer == "16S")
  {
  wine.s <- subset_taxa(wine, !(Domain %in% c("unknown", "Eukaryota",
                                              "Eukaryota_unclassified") |
                                  Phylum %in% c(NA) |
                                  Order %in% c("Chloroplast") |
                                  Family %in% c("Mitochondria")))
  }


unique(tax_table(wine)[, 1])


if(primer == "18S")
  {
  wine.s <- subset_taxa(wine, !(Domain %in% c("Bacteria", "unknown", NA) |
                                  Phylum %in% c("Eukaryota_unclassified",
                                                "Mollusca", "Vertebrata", NA) |
                                  Class %in% c("Insecta", "Ellipura",
                                               "Embryophyta", "Arachnida",
                                               "Heterophyidae", "Ichthyophonae",
                                               "Arthropoda_unclassified",
                                               "unclassified_Hexapoda")))
  }

unique(tax_table(wine.s)[, 2])


################################################################################
### Check for empty taxa and remove if any ####
if(any(taxa_sums(wine.s) == 0))
  {
  sum(taxa_sums(wine.s) == 0)
  wine.s <- prune_taxa(taxa_sums(wine.s) > 0, wine.s)
  }


## Check sample_sums (to check if we should remove a sample with a very low
## number of reads before rarefying)
sums <- sample_sums(wine.s)


barplot(sums, beside = TRUE, col = c("grey"),
        cex.axis = 1, cex.names = 0.6, las = 2)
summary(sums) # looks good, no sample should be removed


## Rarefy
set.seed(100)
wine.r <- rarefy_even_depth(wine.s)
wine.s
wine.r


################################################################################
### Percentage abundance ####
wine.s <- transform_sample_counts(wine.s, function(otu) {otu / sum(otu)})


plot(rowSums(otu_table(wine.s)), ylim = c(0, 1),
     xlab = "Samples", ylab = "Abundance [%]") # 100 % in all samples


################################################################################
### Abundance filtering ####
## For 16S OTU
wine.a <- filter_taxa(wine.s, function(otu) {mean(otu) > 0.000001},
                      prune = TRUE)

## For 16S ASV
# wine.a <- filter_taxa(wine.s, function(otu) {mean(otu) > 0.0001},
#                       prune = TRUE)

points(rowSums(otu_table(wine.a)), col = "red")


if(any(taxa_sums(wine.a) == 0))
  {
  sum(taxa_sums(wine.a) == 0)
  wine.a <- prune_taxa(taxa_sums(wine.a) > 0, wine.a)
  }


wine
wine.s
wine.r
wine.a
otu_table(wine.a)[1:5, 1:5]


################################################################################
### The effect of pruning, rarefying and filtering ####
nr.taxa <- data.frame(NrTaxa = c(dim(tax_table(wine))[1],
                                 dim(tax_table(wine.s))[1],
                                 dim(tax_table(wine.r))[1],
                                 dim(tax_table(wine.a))[1]),
                      Dataset = c("Raw", "Taxonomic filtering", "Rarefaction",
                                  "Abundance filtering"))
nr.taxa$Dataset <- ordered(nr.taxa$Dataset, levels = c("Raw",
                                                       "Taxonomic filtering",
                                                       "Rarefaction",
                                                       "Abundance filtering"))
nr.taxa


## Plot
ggplot(nr.taxa, aes(x = Dataset, y = NrTaxa)) +
  geom_point() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste(alpha, "-diversity"))) +
  xlab("")


################################################################################
## Colorblind scale
cols <- c("#C59434", "#999999", "#009E73") # brown, gray, green


### Plot taxa ####
plot_bar(wine, x = "Phylum", fill = "treatment") +
  facet_grid(Domain ~ ., scales = "free_y", space = "free") +
  geom_bar(aes(color = treatment, fill = treatment), stat = "identity",
           position = "stack") +
  theme(legend.position = "top", legend.direction = "horizontal") +
  xlab("") +
  ylab("Abundance [%]") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  coord_flip()


## If you want to save data tables separately, write raw data tables
wine.a

write.table((otu_table(wine.a)), "wine.a_OTU.csv", col.names = NA,
            row.names = TRUE, sep = ",")
write.table((tax_table(wine.a)), "wine.a_TAX.csv", col.names = NA,
            row.names = TRUE, sep = ",")
write.table((sample_data(wine.a)), "wine.a_meta.csv", col.names = NA,
            row.names = TRUE, sep = ",")


## Extract OTU table with taxonomic annotation
tax <- as.data.frame(wine.a@tax_table@.Data)
otu <- as.data.frame(t(otu_table(wine.a)))

tax.otu <- merge(tax, otu, by = 0, all = TRUE) # by = 0 = by rownames
rownames(tax.otu) <- tax.otu$Row.names
tax.otu$Row.names <- NULL


rm(tax, otu)


### Alpha diversity with untransformed data
## Look at alpha diversity indices
# rarified data
plot_richness(wine.r, measures = c("Observed", "Shannon", "Chao1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste(alpha, "-diversity index"))) +
  xlab("Samples")


# Or grouped by treatment
plot_richness(wine.r, x = "treatment",
              measures = c("Observed", "Shannon", "Chao1"),
              color = "treatment") +
  geom_boxplot(color = "black", alpha = 0.1, outlier.shape = NA) +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text.x = element_blank()) +
  scale_color_manual(values = cols) +
  ylab("Diversity index") +
  xlab("")


## Extract diversity indices as numbers and save as data frame
wine_alpha <- estimate_richness(wine.r, split = TRUE,
                                measure = c("Observed", "Shannon", "Chao1",
                                            "Simpson"))
wine_alpha <- as.data.frame(wine_alpha)
head(wine_alpha)


## Combine metadata with indices data frames for linear regression
reg.df <- cbind(metadata, wine_alpha)
head(reg.df)


### Linear regression ####
library("lme4")


## Look at our data
ggplot(reg.df, aes(x = vineyard, y = Observed, color = treatment)) +
  geom_point(size = 3) +
  ggtitle("Observed") +
  scale_color_manual(values = cols) +
  ylab(expression(alpha-diversity)) +
  xlab("Vineyards")


## Same summarized (caution: use height = 0 to jitter only horizontally!)
ggplot(reg.df, aes(x = treatment, y = Observed, color = treatment)) +
  geom_boxplot(color = "black", outlier.shape = NA) + # omit outliers
  geom_jitter(height = 0, width = 0.3, size = 3) +
  ggtitle("Observed") +
  scale_color_manual(values = cols) +
  theme(legend.position = "none") +
  ylab(expression(alpha-diversity)) +
  xlab("")


### General linear models
## Treatment
m0 <- lmer(Observed ~ 1 + (1|vineyard), data = reg.df, REML = FALSE,
           na.action = "na.fail")

m1 <- lmer(Observed ~ treatment + (1|vineyard) , data = reg.df, REML = FALSE,
           na.action = "na.fail")
summary(m1)


m2 <- lmer(Observed ~ treatment + scale(som) + scale(Cu) + scale(plant_spec) +
             (1|vineyard) , data = reg.df, REML = FALSE, na.action = "na.fail")
summary(m2)


## Model comparison
library("car")


AIC(m0, m1, m2)
anova(m0, m1, m2)


## m2 fits best -> treatment does not explain "Observed" OTU richness. But
## copper does.
m3 <- lmer(Observed ~ scale(Cu) + (1|vineyard) , data = reg.df, REML = FALSE,
           na.action = "na.fail")
summary(m3)
anova(m2, m3)


## Check model fit
res1 <- resid(m3)
qqnorm(res1)
qqline(res1)


## Plot model predictions from m3 Observed ~ Cu
fixef(m3) # look at the fixed effects in m3


ggplot(reg.df, aes(x = scale(Cu), y = Observed, color = treatment)) +
  ggtitle("Observed") +
  geom_point(shape = 16, size = 2) +
  geom_abline(aes(intercept = `(Intercept)`, slope = -67.21558),
              as.data.frame(t(fixef(m3)))) +
  scale_color_manual(values = cols) +
  ylab(expression(alpha-diversity)) +
  xlab("Copper (scaled)") # r: do we have a unit for copper? M:the variable is scaled so guess assigning a unit does not make sense


## Exercise: can you build models for other diversity indices (Chao1, se.chao1,
## Shannon, Simpson)


################################################################################
### Community composition ####
wine.a

## with relative abundance data
otu_table(wine.a)[1:5, 1:5]


## Testing for homogeneous variation of groups in order to exclude significant
## effects due to inhomogeneous distributions
## Convert into distance matrix
d <- distance(wine.a, "bray")


## Homogeneity of dispersion test
sampledf <- data.frame(sample_data(wine.a))

beta <- betadisper(d, sampledf$treatment)
permutest(beta) # not significant: variations are homogeneous


## Perform test
anova(beta)


## Visualize variances
beta <- with(sampledf, betadisper(d, treatment))
plot(beta)
boxplot(beta, xlab = "Treatment", col = cols)


################################################################################
### Unconstrained ordination ####
wine.log <- transform_sample_counts(wine.a, function(otu) {log1p(otu)})



## NMDS of Bray-Curtis distance
wine.nmds <- ordinate(wine.log, "NMDS", "bray", autotransform = FALSE,
                      # weakties = FALSE, # activate if stress is (nearly) zero
                      trymax = 50)
wine.nmds

stressplot(wine.nmds) # GOF


plot_ordination(wine.log, wine.nmds, color = "treatment", shape = "treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  scale_color_manual(values = cols) +
  geom_text(aes(label = sample_data(wine.a)$vineyard), color = "black",
            size = 3) +
  ggtitle("NMDS of Bray-Curtis distance") +
  coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = treatment), type = "t", linetype = 2, size = 0.2)


## Remove outlier sample (if justified)
wine.log.out <- subset_samples(wine.log, vineyard != "31") #2015 = "x139-15"


## NMDS of Bray-Curtis distance without outliers
wine.nmds <- ordinate(wine.log.out, "NMDS", "bray")
wine.nmds

stressplot(wine.nmds)


plot_ordination(wine.log.out, wine.nmds, color = "treatment",
                shape = "treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  scale_color_manual(values = cols) +
  geom_text(aes(label = sample_data(wine.log.out)$vineyard), color = "black",
            size = 3) +
  ggtitle("NMDS of Bray-Curtis distance outliers removed") +
  coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = treatment), type = "t", linetype = 2, size = 0.2)


### MDS ###
wine.mds <- ordinate(wine.log, method = "PCoA", distance = "bray")


barplot(wine.mds$values$Relative_eig)
biplot(wine.mds, data.frame(otu_table(wine.log)))


plot_ordination(wine.log, wine.mds, shape = "treatment",
                color = "treatment", title = NULL, label = "vineyard",
                axes = 1:2) +
  stat_ellipse(aes(group = treatment), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = cols) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box="vertical")


################################################################################
### Constrained ordination ####
## RDA
wine.rda <- ordinate(wine.log, "RDA", "bray")

ordcap <- ordinate(wine.log, "CAP", "bray", ~ treatment + Cu + som)
summary(ordcap)


plot_ordination(wine.log, ordcap, "samples", color = "treatment") +
  scale_color_manual(values = cols) +
  coord_fixed(ratio = 1)


# TO DO #### RDA with environmental variables (Cu, som, plant diversity)
