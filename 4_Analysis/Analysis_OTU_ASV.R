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
library("microbiome")


rm(list = ls())
setwd("/home/rstudio/")


################################################################################
set.seed(11948)


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
method <- "OTU"
primer <- "16S"
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


MetaData <- sample_data(metadata)
wine <- merge_phyloseq(wine, MetaData)


## Exploring our data
nsamples(wine)
ntaxa(wine)
sample_names(wine)[1:66]
sample_variables(wine) # metadata variables
sample_data(wine)
otu_table(wine)[1:5, 1:5]
tax_table(wine)[1:5, 1:6]


# OR
summarize_phyloseq(wine)


## Which taxa do we have in our data? Here on the phyla level:
get_taxa_unique(wine, taxonomic.rank = rank_names(wine)[1], errorIfNULL = TRUE)


## Look at our data
readsumsdf <- data.frame(nreads = taxa_sums(wine),
                         sorted = 1:ntaxa(wine), type = "OTUs")


## Plot our raw reads per OTU
ggplot(readsumsdf, aes(x = sorted, y = nreads)) +
  geom_bar(stat = "identity") +
  ggtitle("Total number of reads") +
  scale_y_log10() +
  # coord_cartesian(ylim = c(0, 100)) +
  xlab("") +
  ylab(expression(paste(log[10], "(Number of reads)")))


wine
colSums(readsumsdf[, 1, drop = FALSE])

################################################################################
### Remove singletons ###
sum(taxa_sums(wine) == 1)


if(any(taxa_sums(wine) == 1))
  {wine <- prune_taxa(taxa_sums(wine) > 1, wine)}
wine


################################################################################
### Remove spurious taxa ####
if(primer == "16S")
  {
  wine.s <- subset_taxa(wine, !(Domain %in% c("unknown", "Eukaryota") |
                                  Phylum %in% c("Eukaryota_unclassified",
                                                NA) |
                                  Order %in% c("Chloroplast") |
                                  Family %in% c("Mitochondria")))
  }


if(primer == "18S")
  {
  wine.s <- subset_taxa(wine, !(Domain %in% c("Bacteria", "unknown") |
                                  Phylum %in% c("Eukaryota_unclassified",
                                                "Mollusca", "Vertebrata",
                                                NA) |
                                  Class %in% c("Insecta", "Ellipura",
                                               "Embryophyta", "Arachnida",
                                               "Heterophyidae",
                                               "Ichthyophonae",
                                               "Arthropoda_unclassified",
                                               "unclassified_Hexapoda")))
  }


wine
wine.s


################################################################################
### Check for empty taxa and remove if any ####
sum(taxa_sums(wine.s) == 0)

if(any(taxa_sums(wine.s) == 0))
  {wine.s <- prune_taxa(taxa_sums(wine.s) > 0, wine.s)}


## Check sample_sums (to check if we should remove a sample with a very low
## number of reads before rarefying)
sums <- sample_sums(wine.s)


barplot(sums, beside = TRUE, col = c("grey"),
        cex.axis = 1, cex.names = 0.6, las = 2)
summary(sums) # looks good no sample should be removed


## Rarefy
# set.seed(100)
wine.r <- rarefy_even_depth(wine)
wine.r
wine.s


################################################################################
### Percentage abundance ####
wine.s <- transform_sample_counts(wine.s, function(otu) {otu / sum(otu)})
plot(rowSums(otu_table(wine.s)), ylim = c(0, 1),
     xlab = "Samples", ylab = "Abundance [%]") # 100 % in all samples


################################################################################
### Abundance filtering ####
wine.a <- filter_taxa(wine.s, function(otu) {mean(otu) > 0.00001},
                      prune = TRUE)
points(rowSums(otu_table(wine.a)), col = "red")


sum(taxa_sums(wine) == 0)
if(any(taxa_sums(wine) == 0))
  {wine <- prune_taxa(taxa_sums(wine) > 0, wine)}


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
## Colorblind scale
cols <- c("#C59434", "#999999", "#009E73") # brown, gray, green


### Plot taxa ####
plot_bar(wine.a, x = "Phylum", fill = "ord.treatment") +
  facet_grid(Domain ~ ., scales = "free_y", space = "free") +
  geom_bar(aes(color = ord.treatment, fill = ord.treatment), stat = "identity",
           position = "stack") +
  theme(legend.position = "top", legend.direction = "horizontal") +
  xlab("") +
  ylab("Abundance [%]") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  coord_flip()


### Plot rank abundance of Orders ####
colnames(tax_table(wine.a)) # print the available taxonomic ranks
wine.ord <- aggregate_taxa(wine.a, "Order")
wine.ord


par(mar = c(10, 4, 4, 2) + 0.1)  # make more room on bottom margin
N <- 20 # show the top20 abundant Orders
barplot(sort(taxa_sums(wine.ord), TRUE)[1:N] / nsamples(wine.ord),
        main = "Relative abundance of top 20 most abundant Orders",
        las = 2)


## If you want to save data tables separately, write raw data tables
wine.a

write.table((otu_table(wine.a)), "wine.a_OTU.csv", col.names = NA,
            row.names = TRUE, sep = ",")
write.table((tax_table(wine.a)), "wine.a_TAX.csv", col.names = NA,
            row.names = TRUE, sep = ",")
write.table((sample_data(wine.a)), "wine.a_meta.csv", col.names = NA,
            row.names = TRUE, sep = ",")


## Extract OTU table with taxonomic annotation into a data.frame
tax <- as.data.frame(wine.a@tax_table@.Data)
otu <- as.data.frame(t(otu_table(wine.a)))

tax.otu <- merge(tax, otu, by = 0, all = TRUE) # by = 0 = by rownames
rownames(tax.otu) <- tax.otu$Row.names
tax.otu$Row.names <- NULL
tax.otu[1:5, 1:10]

rm(tax, otu)


### Alpha diversity with untransformed data
## Look at alpha diversity indices
## Rarified data

plot_richness(wine.r, measures = c("Observed", "Shannon", "Chao1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste(alpha, "-diversity index"))) +
  xlab("Samples")


# Or by treatment
plot_richness(wine.r, x = "ord.treatment",
              measures = c("Observed", "Shannon", "Chao1"),
              color = "treatment") +
  geom_boxplot(color = "black", alpha = 0.1, outlier.shape = NA) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
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
ggplot(reg.df, aes(x = vineyard, y = Observed, color = ord.treatment)) +
  geom_point(size = 3) +
  ggtitle("Observed") +
  scale_color_manual(values = cols) +
  ylab(expression(alpha-diversity)) +
  xlab("Vineyards")


## Same summarized (caution: use height = 0 to jitter only horizontally!)
ggplot(reg.df, aes(x = ord.treatment, y = Observed, color = treatment)) +
  geom_boxplot(color = "black", outlier.shape = NA) + # omit outliers
  geom_jitter(height = 0, width = 0.3, size = 3) +
  ggtitle("Observed") +
  scale_color_manual(values = cols) +
  theme(legend.position = "none") +
  ylab(expression(alpha-diversity)) +
  xlab("")


## Set treatment variable as.factor
reg.df$treatment <- as.factor(reg.df$treatment)
str(reg.df)


## Relevel default reference for comparison "BareGround"
reg.df <- within(reg.df, treatment <- relevel(treatment, ref = "BareGround"))
str(reg.df)


### General linear models
## Treatment
m0 <- lmer(Observed ~ 1 + (1|vineyard), data = reg.df, REML = FALSE,
           na.action = "na.fail")

m1 <- lmer(Observed ~ treatment + (1|vineyard) , data = reg.df, REML = FALSE,
           na.action = "na.fail")
summary(m1) # t value is higher (+ or -) or equal to 1.96 it is considered significant (rule-of-thumb)


m2 <- lmer(Observed ~ treatment + scale(Cu) +
             (1|vineyard) , data = reg.df, REML = FALSE, na.action = "na.fail")
summary(m2)


m3 <- lmer(Observed ~ treatment + scale(som) + scale(Cu) + scale(plant_spec) +
             (1|vineyard) , data = reg.df, REML = FALSE, na.action = "na.fail")
summary(m3)


## Model comparison
library("car")


AIC(m0, m1, m2, m3)
anova(m0, m1, m2, m3)

# The best model seems to be m2.
# m2 fits best -> treatment does not explain "Observed" OTU richness alone. Copper shows a larger negative effect on "Observed".


## Check model fit
res1 <- resid(m2)
qqnorm(res1)
qqline(res1)
plot(res1)


## Plot model predictions from m2 Observed ~ Cu
fixef(m2) # look at the fixed effects in m2
slope <- fixef(m2)[4] # extract slope of Cu for regression plot

ggplot(reg.df, aes(x = scale(Cu), y = Observed, color = ord.treatment)) +
  ggtitle("Observed") +
  geom_point(shape = 16, size = 2) +
  geom_abline(aes(intercept = `(Intercept)`, slope = slope),
              as.data.frame(t(fixef(m2)))) +
  scale_color_manual(values = cols) +
  ylab(expression(alpha-diversity)) +
  xlab("Copper (scaled)") # Copper is in ??g/kg soil but here the variable is scaled in the analysis so assigning a unit does not make sense


## Exercise: can you build models for other diversity indices (Chao1, se.chao1,
## Shannon, Simpson)


################################################################################
### Community composition ####
wine.a


## With relative abundance data
otu_table(wine.a)[1:5, 1:5]


## Testing for homogeneous variation of groups in order to exclude significant
## effects due to inhomogeneous distributions
## Convert into distance matrix
d <- distance(wine.a, "bray")


## Homogeneity of dispersion test
sampledf <- data.frame(sample_data(wine.a))

beta <- betadisper(d, sampledf$treatment)
permutest(beta) # not significant: variations are sufficiently homogeneous


## Perform test
anova(beta)


## Visualize variances
beta <- with(sampledf, betadisper(d, treatment))
plot(beta)
boxplot(beta, xlab = "Treatment", col = cols)


### NMDS (unconstrained Ordination) ####
## NMDS of Bray-Curtis distance
p_nmds <- ordinate(wine.a, "NMDS", "bray", autotransform = FALSE,
                   trymax = 50)
p_nmds

stressplot(p_nmds) # GOF


plot_ordination(wine.a, p_nmds, color = "ord.treatment", shape = "ord.treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  scale_color_manual(values = cols) +
  geom_text(aes(label = sample_data(wine.a)$vineyard), color = "black",
            size = 3) +
  ggtitle("NMDS of Bray-Curtis distance") +
  coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = ord.treatment), type = "t", linetype = 2, size = 0.2)


## Remove outlier sample (if justified)
wine.a.out <- subset_samples(wine.a, vineyard != "31")


## NMDS of Bray-Curtis distance outliers removed
p_nmds <- ordinate(wine.a.out, "NMDS", "bray")
p_nmds

stressplot(p_nmds)


plot_ordination(wine.a.out, p_nmds, color = "ord.treatment",
                shape = "ord.treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  scale_color_manual(values = cols) +
  geom_text(aes(label = sample_data(wine.a.out)$vineyard), color = "black",
            size = 3) +
  ggtitle("NMDS of Bray-Curtis distance outliers removed") +
  coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = ord.treatment), type = "t", linetype = 2, size = 0.2)


### Constrained ordination ####
## RDA

## RDA
# p_rda <- ordinate(wine.a, "RDA", "bray")
#
# ordcap <- ordinate(wine.a, "CAP", "bray", ~ treatment + Cu + som)
# summary(ordcap)
#
#
# plot_ordination(wine.a, ordcap, "samples", color = "treatment") +
#   scale_color_manual(values = cols)

## CAP (dbRDA) ordinate ####
cap_ord <- ordinate(
  physeq = wine.a,
  method = "CAP",
  distance = "bray",
  formula = ~ treatment + Cu + som + plant_spec)
summary(cap_ord)
cap_ord
anova(cap_ord)
screeplot(cap_ord) # visualisation of explained variation by axes

RsquareAdj(cap_ord)


## CAP plot ####
cap_plot <- plot_ordination(
  physeq = wine.a,
  ordination = cap_ord,
  color = "treatment",
  shape = "treatment",
  axes = c(1,2))

cap_plot


## Customize your plot
cap_plot = cap_plot + geom_point(size = 4) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c(18, 16, 17)) +
  ggtitle("Community composition dbRDA")
#  geom_text(aes(label=sample_data(wine.a)$vineyard), color = "black", size = 2.5) #enable if you wnt your sites to be numbered

cap_plot

## Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

## Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

## Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1,
                 yend = CAP2,
                 x = 0,
                 y = 0,
                 shape = NULL,
                 color = NULL,
                 label = labels)

label_map <- aes(x = 1.3 * CAP1,
                 y = 1.3 * CAP2,
                 shape = NULL,
                 color = NULL,
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))


# Make a new graphic
p <- cap_plot +
  geom_segment(
    mapping = arrow_map,
    size = .8,
    data = arrowdf,
    color = "grey28",
    arrow = arrowhead
  ) +
  geom_text(
    mapping = label_map,
    size = 4,
    data = arrowdf,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")
p


### The core microbiome - common core across all treatments ####
wine.core <- core(wine.a, detection = .0005, prevalence = .90) # a minimum abundance of 0.05 % and prevalent in 90 % of the samples
summarize_phyloseq(wine.core)
summarize_phyloseq(wine.a)
wine.core
wine.a


## Percentage of core n_taxa compared to full dataset
ntaxa(wine.core) / ntaxa(wine.a) * 100


### Plot Venn Diagrams ###
library("VennDiagram")


## Prepare datasets for VennDiagramm
## Bare ground
all.BG <- subset_samples(wine.a, treatment == "BareGround") # select bare ground samples
any(taxa_sums(all.BG) == 0) # OTUs with abundance 0
all.BG <- prune_taxa(taxa_sums(all.BG) > 0, all.BG) # only keep OTUs >0
wine.BG <- row.names(tax_table(all.BG))


## Alternating cover
all.AC <- subset_samples(wine.a, treatment == "AlternatingCover") # select bare ground samples
any(taxa_sums(all.AC) == 0) # OTUs with abundance 0
all.AC <- prune_taxa(taxa_sums(all.AC) > 0, all.AC) # only keep OTUs >0
wine.AC <- row.names(tax_table(all.AC))


## Complete cover
all.CC <- subset_samples(wine.a, treatment == "CompleteCover") # select bare ground samples
any(taxa_sums(all.CC) == 0) # OTUs with abundance 0
all.CC <- prune_taxa(taxa_sums(all.CC) > 0, all.CC) # only keep OTUs >0
wine.CC <- row.names(tax_table(all.CC))


all <- list(wine.BG, wine.AC, wine.CC)


## VennDiagram for all OTUs
dev.off()
plt <- venn.diagram(
  all,
  category.names = c("Bare ground" , "Alternating cover" , "Complete cover"),
  filename = NULL)
grid::grid.draw(plt)


## Prepare datasets for VennDiagram with core microbiome per treatment
## Bare ground
BG <- subset_samples(wine.a, treatment == "BareGround") # select bare ground samples
core.BG <- wine.core <- core(BG, detection = .0005, prevalence = .90) # a minimum abundance of 0.05 % and prevalent in 95% of the samples
core.BG <- prune_taxa(taxa_sums(core.BG) > 0, core.BG) # only keep OTUs >0
core.BG <- row.names(tax_table(core.BG))


## Alternating cover
AC <- subset_samples(wine.a, treatment == "AlternatingCover") # select bare ground samples
core.AC <- wine.core <- core(AC, detection = .0005, prevalence = .90) # a minimum abundance of 0.05 % and prevalent in 95% of the samples
core.AC <- prune_taxa(taxa_sums(core.AC) > 0, core.AC) # only keep OTUs >0
core.AC <- row.names(tax_table(core.AC))


## Complete cover
CC <- subset_samples(wine.a, treatment == "CompleteCover") # select bare ground samples
core.CC <- wine.core <- core(CC, detection = .0005, prevalence = .90) # a minimum abundance of 0.05 % and prevalent in 95% of the samples
core.CC <- prune_taxa(taxa_sums(core.CC) > 0, core.CC) # only keep OTUs >0
core.CC <- row.names(tax_table(core.CC))

core <- list(core.BG, core.AC, core.CC)


## Chart
dev.off()
plt <- venn.diagram(
  core,
  category.names = c("Bare ground" , "Alternating cover" , "Complete cover"),
  filename = NULL)
grid::grid.draw(plt)


####
dev.off()


plt <- venn.diagram(core,
                    category.names = c("Bare ground" , "Alternating cover" ,
                                       "Complete cover"),
                    filename = NULL,

                    # Circles
                    lwd = 2,
                    lty = 'blank',
                    fill = cols,

                    # Numbers
                    cex = .6,
                    fontface = "bold",
                    fontfamily = "sans",

                    # Set names
                    cat.cex = 0.6,
                    cat.fontface = "bold",
                    cat.default.pos = "outer",
                    cat.pos = c(-27, 27, 135),
                    cat.dist = c(0.055, 0.055, 0.085),
                    cat.fontfamily = "sans",
                    rotation = 1)
grid::grid.draw(plt)


################################################################################
################################################################################



