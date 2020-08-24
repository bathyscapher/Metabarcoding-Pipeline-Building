################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### mothur pipeline OTU: rarefaction curves
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2020
### korn@cumulonimbus.at
################################################################################
################################################################################
# set working directory ####

library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("vegan")


rm(list = ls())


setwd("C:/Users/Steima/Seafile/My Library/Doctoral school CUSO/Cuso organising a course/2020/Workshop-Analysis/")
setwd("~/Metabarcoding-Pipeline-Building/")


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
### Read taxa into phyloseq object ####
## Note: for now choose RDP for OTU, DECIPHER for ASV
method <- "OTU"
primer <- "16S"
classifier <- "RDP"
# classifier <- "DECIPHER"


wine <- readTaxa(method = method, primer = primer, classifier = classifier)
wine


colnames(tax_table(wine))


################################################################################
### Add metadata ####
## Load metadata
metadata <- read.csv("Metadata_corr.csv", sep = ",", header = TRUE,
                     row.names = 1)
head(metadata)


## Rename and sort factors
metadata$treatment <- factor(metadata$treatment)
levels(metadata$treatment) <- c("Alternating cover", "Bare ground",
                                "Complete cover")
metadata$treatment <- ordered(metadata$treatment,
                              levels = c("Bare ground", "Alternating cover",
                                         "Complete cover"))


MetaData <- sample_data(metadata)
wine <- merge_phyloseq(wine, MetaData)


#M
## Exploring our data
nsamples(wine)
ntaxa(wine)
sample_names(wine)[1:66] # why the subset if it includes the full dataset?

sample_variables(wine) # metadata variables
sample_data(wine)
otu_table(wine)[1:5, 1:5]
tax_table(wine)[1:5, 1:6]


## Which taxa do we have in our data?
get_taxa_unique(wine, taxonomic.rank = rank_names(wine)[2], errorIfNULL = TRUE)

## Look at our data
readsumsdf <- data.frame(nreads = sort(taxa_sums(wine), TRUE),
                         sorted = 1:ntaxa(wine), type = "OTUs")

## Plot our raw reads per OTU
p <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) +
  geom_bar(stat = "identity") +
  ggtitle("Total number of reads") +
  scale_y_log10() +
  facet_wrap(~ type, 1, scales = "free") + # scales = "free" has no function here, remove?
  xlab("") +
  ylab(expression(paste(log[10], "(Number of reads)")))
p


wine
colSums(readsumsdf[, 1, drop = FALSE])


################################################################################
### Remove spurious taxa ####
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
### Check for empty taxa and remove if any ####
if(any(taxa_sums(wine) == 0))
{
  sum(taxa_sums(wine) == 0)
  wine <- prune_taxa(taxa_sums(wine) > 0, wine)
}

#M
## Check sample_sums (to check if we should remove a sample with a very low
## number of reads before rarefying)
sums <- sample_sums(wine)

barplot(sums, beside = TRUE, col = c("grey"),
        cex.axis = 1, cex.names = 0.6, las = 2)
summary(sums) # looks good no sample should be removed


## Rarefy
set.seed(100)
wine.r <- rarefy_even_depth(wine.s)
wine.r #
wine


################################################################################
### Percentual abundance ####
wine.s <- transform_sample_counts(wine.r, function(otu) {otu / sum(otu)})
plot(rowSums(otu_table(wine.s)), ylim = c(0, 1),
     xlab = "Samples", ylab = "Abundance [%]") # 100 % in all samples


################################################################################
### Abundance filtering ####
wine.a <- filter_taxa(wine.s, function(otu) {mean(otu) > 0.00001},
                      prune = TRUE)
points(rowSums(otu_table(wine.a)), col = "red")


if(any(taxa_sums(wine) == 0))
{
  sum(taxa_sums(wine) == 0)
  wine <- prune_taxa(taxa_sums(wine) > 0, wine)
}

wine.s
wine.a
otu_table(wine.a)[1:5, 1:5]


################################################################################
### The effect of pruning, rarefying and filtering ####
nr.taxa <- data.frame(NrTaxa = c(dim(tax_table(wine))[1],
                                 dim(tax_table(wine.r))[1],
                                 dim(tax_table(wine.s))[1],
                                 dim(tax_table(wine.a))[1]),
                      Dataset = c("Raw", "Rarefaction", "Taxonomic filtering",
                                  "Abundance filtering"))
nr.taxa$Dataset <- ordered(nr.taxa$Dataset, levels = c("Raw", "Rarefaction",
                                                       "Taxonomic filtering",
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
### Plot taxa ####
plot_bar(wine.a, x = "Phylum", fill = "Phylum") +
  facet_grid(Domain ~ ., scales = "free_y", space = "free") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity",
           position = "stack") +
  theme(legend.position = "top", legend.direction = "horizontal") +
  xlab("") +
  ylab("Abundance [%]") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  coord_flip()


#M
#If you want to save data tables separately
## Write rawdata tables
wine.a

write.table((otu_table(wine.a)), "wine.a_OTU.csv", col.names = NA,
            row.names = TRUE, sep = ",")
write.table((tax_table(wine.a)), "wine.a_TAX.csv", col.names = NA,
            row.names = TRUE, sep = ",")
write.table((sample_data(wine.a)), "wine.a_meta.csv", col.names = NA,
            row.names = TRUE, sep = ",")


### Alpha diversity with untransformed data
## Look at alpha diversity indices
# raw data. r: that's the rarefied data, right?
plot_richness(wine.r, measures = c("Observed", "Shannon", "Chao1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste(alpha, "-diversity index"))) +
  xlab("Samples")


# r: Or like this? (Chao1 would need a jitter if we keep it...)
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

## Colorblind scale
cols <- c("#C59434", "#999999", "#009E73") # brown, gray, green


## Look at my data
Obs <- ggplot(reg.df, aes(x = vineyard, y = Observed, color = treatment)) +
  geom_point(size = 3) +
  ggtitle("Observed") +
  scale_color_manual(values = cols) +
  ylab(expression(alpha-diversity)) +
  xlab("Vineyards")
Obs


## Same summarized (caution: use height = 0 to jitter only horizontally!)
ggplot(reg.df, aes(x = treatment, y = Observed, color = treatment)) +
  geom_boxplot(color = "black", outlier.shape = NA) + # omit outliers
  geom_jitter(height = 0, size = 3) +
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


## m2 fits best -> treatment does not expain "Observed" OTU richness. But copper
## does.
m3 <- lmer(Observed ~ scale(Cu) + (1|vineyard) , data = reg.df, REML = FALSE,
           na.action = "na.fail")
summary(m3)
anova(m2, m3)


## Check modelfit
res1 <- resid(m3)
qqnorm(res1)
qqline(res1)


## Plot model predictions from m3 Observed ~ Cu
fixef(m3) #look at the fixed effects in m3


obs <- ggplot(reg.df, aes(x = scale(Cu), y = Observed, color = treatment)) +
  ggtitle("Observed") +
  geom_point(shape = 16, size = 2) +
  geom_abline(aes(intercept = `(Intercept)`, slope = -67.21558),
              as.data.frame(t(fixef(m3)))) +
  scale_color_manual(values = cols) +
  ylab(expression(alpha-diversity)) +
  xlab("Copper [unit]") # do we have a unit for copper?
obs


# You can build models for other diversity indices (Chao1, se.chao1, Shannon, Simpson)


################################################################################
### Community composition ####
wine.a

## with relative abundance data
otu_table(wine.a)[1:5, 1:5]

## Testing for homogenous variation of groups in order to exclude significant
## effect due to inhomogenous ditributions
## Convert into distance matrix
d <- distance(wine.a, "bray")


## Homogeneity of dispersion test
sampledf <- data.frame(sample_data(wine.a))

beta <- betadisper(d, sampledf$treatment)
permutest(beta) # not significant: variations are homogenous


## Perform test
anova(beta)


## visualize variances
beta <- with(sampledf, betadisper(d, treatment))
plot(beta)
boxplot(beta, xlab = "Treatment")


### NMDS (unconstrained Ordination) ####
## NMDS of Bray-Curtis distance
# r: maybe better set autotransform to FALSE and select a transformation?
p_nmds <- ordinate(wine.a, "NMDS", "bray", autotransform = TRUE)
p_nmds

stressplot(p_nmds)


p <- plot_ordination(wine.a, p_nmds, color = "treatment", shape = "treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  # scale_color_manual(values = c('#CCCC00', '#996600','#66CC33')) +
  scale_color_manual(values = cols) +
  geom_text(aes(label=sample_data(wine.a)$vineyard), color = "black", size = 3) +
  ggtitle("NMDS of Bray-Curtis distance") +
  coord_fixed(ratio = 1)
p
p + stat_ellipse(aes(group = treatment), type = "t", linetype = 2, size = 0.2)


## Remove outlier sample (if justified)
wine.a.out <- subset_samples(wine.a, vineyard != "31") #2015 = "x139-15"


## NMDS of Bray-Curtis distance outliers removed
p_nmds <- ordinate(wine.a.out, "NMDS", "bray")
p_nmds

stressplot(p_nmds)


p <- plot_ordination(wine.a.out, p_nmds, color = "treatment",
                     shape = "treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  # scale_color_manual(values = c('#CCCC00', '#996600','#66CC33')) +
  scale_color_manual(values = cols) +
  geom_text(aes(label = sample_data(wine.a.out)$vineyard), color = "black",
            size = 3) +
  ggtitle("NMDS of Bray-Curtis distance outliers removed") +
  coord_fixed(ratio = 1)
p
p + stat_ellipse(aes(group = treatment), type = "t", linetype = 2, size = 0.2)


### Constrained ordination ####
## RDA
p_rda <- ordinate(wine.a.out, "RDA", "bray")

ordcap <- ordinate(wine.a, "CAP", "bray", ~ treatment + Cu + som)
summary(ordcap)


plot_ordination(wine.a, ordcap, "samples", color = "treatment") +
  scale_color_manual(values = cols)


# TO DO #### RDA with environmental variables (Cu, som, plant diversity)


### MDS = PCoA ####
## First, log-transform
wine.log <- transform_sample_counts(wine.a, function(otu) {log1p(otu)})

wine.nmds <- ordinate(wine.log, method = "PCoA", distance = "bray",
                      autotransform = FALSE)

barplot(wine.nmds$values$Relative_eig)
biplot(wine.nmds, data.frame(otu_table(wine.a)))


plot_ordination(wine.a, wine.nmds, color = "treatment", title = NULL) +
  stat_ellipse(aes(group = treatment), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  scale_color_manual(values = cols) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal")


################################################################################
################################################################################
