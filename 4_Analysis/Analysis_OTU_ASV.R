###############################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### Analysis of microbiome data
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2020
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


setwd("C:/Users/Steima/Seafile/My Library/Doctoral school CUSO/Cuso organising a course/2020/Workshop-Analysis/") # set working directory to home


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

## Which taxa do we have in our data? Here on the phlya level:
get_taxa_unique(wine, taxonomic.rank = rank_names(wine)[2], errorIfNULL = TRUE)


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
### Remove spurious taxa ####
if(primer == "16S")
{
  wine.s <- subset_taxa(wine, !(Domain %in% c("unknown", "Eukaryota") |
                                  Phylum %in% c("Eukaryota_unclassified", NA) |
                                  Order %in% c("Chloroplast") |
                                  Family %in% c("Mitochondria")))
}

if(primer == "16S")
{
  wine.s <- subset_taxa(wine, !(Domain %in% c("Bacteria", "unknown") |
                                  Phylum %in% c("Eukaryota_unclassified",
                                                "Mollusca", "Vertebrata", NA) |
                                  Class %in% c("Insecta", "Ellipura",
                                               "Embryophyta", "Arachnida",
                                               "Heterophyidae", "Ichthyophonae",
                                               "Arthropoda_unclassified",
                                               "unclassified_Hexapoda")))
}


################################################################################
### Check for empty taxa and remove if any ####
if(any(taxa_sums(wine.s) == 0))
{
  sum(taxa_sums(wine.s) == 0)
  wine <- prune_taxa(taxa_sums(wine.s) > 0, wine.s)
}


## Check sample_sums (to check if we should remove a sample with a very low
## number of reads before rarefying)
sums <- sample_sums(wine.s)


barplot(sums, beside = TRUE, col = c("grey"),
        cex.axis = 1, cex.names = 0.6, las = 2)
summary(sums) # looks good no sample should be removed


## Rarefy
set.seed(100)
wine.r <- rarefy_even_depth(wine)
wine.r
wine


################################################################################
### Percentual abundance ####
wine.s <- transform_sample_counts(wine.s, function(otu) {otu / sum(otu)})
plot(rowSums(otu_table(wine.s)), ylim = c(0, 1),
     xlab = "Samples", ylab = "Abundance [%]") # 100 % in all samples


################################################################################
### Abundance filtering ####
wine.a <- filter_taxa(wine.s, function(otu) {mean(otu) > 0.0001},
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
## Colorblind scale
cols <- c("#C59434", "#999999", "#009E73") # brown, gray, green


### Plot taxa ####
plot_bar(wine.a, x = "Phylum", fill = "treatment") +
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
  

# plot rank abundance of Orders ####

#library(devtools) # Load the devtools package
#install_github("microbiome/microbiome") # Install the package

colnames(tax_table(wine.a)) # print the available taxonomic ranks
wine.ord <- aggregate_taxa(wine.a, "Order")

par(mar = c(10, 4, 4, 2) + 0.1)  # make more room on bottom margin
N <- 20 #show the top20 abundant Orders
barplot(sort(taxa_sums(wine.ord), TRUE)[1:N]/nsamples(wine.ord), 
        main="Relative abundance of top 20 most abundant Orders",
        las = 2)


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
head(tax.otu)

rm(tax, otu)


### Alpha diversity with untransformed data
## Look at alpha diversity indices
# rarified data
plot_richness(wine.r, measures = c("Observed", "Shannon", "Chao1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste(alpha, "-diversity index"))) +
  xlab("Samples")


# r: Or like this?
plot_richness(wine.r, x = "treatment",
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


## m2 fits best -> treatment does not explain "Observed" OTU richness alone. But
## copper shows a larger negative effect on "Observed".
m3 <- lmer(Observed ~ scale(Cu) + (1|vineyard) , data = reg.df, REML = FALSE,
           na.action = "na.fail")
summary(m3)
anova(m1, m2, m3)
# however the best model ist  still the model that cotains treatment and other variables

## Check model fit
res1 <- resid(m3)
qqnorm(res1)
qqline(res1)


## Plot model predictions from m3 Observed ~ Cu
fixef(m2) # look at the fixed effects in m3


ggplot(reg.df, aes(x = scale(Cu), y = Observed, color = treatment)) +
  ggtitle("Observed") +
  geom_point(shape = 16, size = 2) +
  geom_abline(aes(intercept = `(Intercept)`, slope = -67.21558),
              as.data.frame(t(fixef(m3)))) +
  scale_color_manual(values = cols) +
  ylab(expression(alpha-diversity)) +
  xlab("Copper (scaled)") #  Copper is in µg/kg soil but here the variable is scaled in the analysis so assigning a unit does not make sense


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


plot_ordination(wine.a, p_nmds, color = "treatment", shape = "treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  scale_color_manual(values = cols) +
  geom_text(aes(label = sample_data(wine.a)$vineyard), color = "black",
            size = 3) +
  ggtitle("NMDS of Bray-Curtis distance") +
  coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = treatment), type = "t", linetype = 2, size = 0.2)


## Remove outlier sample (if justified)
wine.a.out <- subset_samples(wine.a, vineyard != "31") 


## NMDS of Bray-Curtis distance outliers removed
p_nmds <- ordinate(wine.a.out, "NMDS", "bray")
p_nmds

stressplot(p_nmds)


plot_ordination(wine.a.out, p_nmds, color = "treatment",
                shape = "treatment") +
  geom_point(size = 5) +
  scale_shape_manual(values = c(18, 16, 17)) +
  scale_color_manual(values = cols) +
  geom_text(aes(label = sample_data(wine.a.out)$vineyard), color = "black",
            size = 3) +
  ggtitle("NMDS of Bray-Curtis distance outliers removed") +
  coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = treatment), type = "t", linetype = 2, size = 0.2)


### Constrained ordination ####
## RDA

# CAP (dbRDA) ordinate ####
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


# CAP plot ####
cap_plot <- plot_ordination(
  physeq = wine.a,
  ordination = cap_ord, 
  color = "treatment", 
  shape = "treatment",
  axes = c(1,2))
cap_plot = cap_plot + geom_point(size = 4) + ggtitle() +
  scale_shape_manual(values=c(18, 16, 17))+
  scale_color_manual(values=c('#996600','#CCCC00', '#66CC33', "#667C33"))
#  geom_text(aes(label=sample_data(ps)$vineyard), color = "black", size = 2.5)

cap_plot 


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
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
    size = .7, 
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
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_bw()

p

#  The core microbiome - common core across all treatments ####

wine.core <- core(wine.a, detection = .0005, prevalence = .90) # a minimum abundance of 0.05 % and prevalent in 95% of the samples
summarize_phyloseq(wine.core)
summarize_phyloseq(wine.a)


wine.core_nreads <- data.frame(nreads = sort(taxa_sums(wine.core), TRUE),
                         sorted = 1:ntaxa(wine.core), type = "OTUs")
colSums(wine.core_nreads[, 1, drop = FALSE])


# percentage of core n_taxa compared to full dataset
ntaxa(wine.core)/ntaxa(wine.a)*100


# plot Venn Diagrams
library("VennDiagram")

# prepare datasets for VennDiagramm
#Bare ground
all.BG <- subset_samples(wine.a, treatment == "Bare ground") # select bare ground samples
any(taxa_sums(all.BG) == 0) # OTUs with abundance 0
all.BG = prune_taxa(taxa_sums(all.BG) > 0, all.BG) # only keep OTUs lager >0
wine.BG <- row.names(tax_table(all.BG))

# Alternating cover 
all.AC <- subset_samples(wine.a, treatment == "Alternating cover") # select bare ground samples
any(taxa_sums(all.AC) == 0) # OTUs with abundance 0
all.AC = prune_taxa(taxa_sums(all.AC) > 0, all.AC) # only keep OTUs lager >0
wine.AC <- row.names(tax_table(all.AC))

# Complete cover 
all.CC <- subset_samples(wine.a, treatment == "Complete cover") # select bare ground samples
any(taxa_sums(all.CC) == 0) # OTUs with abundance 0
all.CC = prune_taxa(taxa_sums(all.CC) > 0, all.CC) # only keep OTUs lager >0
wine.CC <- row.names(tax_table(all.CC))

all=list(wine.BG, wine.AC, wine.CC)

# VennDiagram for all OTUs
venn.diagram(
  all,
  category.names = c("Bare ground" , "Alternating cover" , "Complete cover"),
  filename = 'Venn_all.png',
  output=TRUE
)

# prepare datasets for VennDiagram with core microbiome per treatment 
#Bare ground
BG <- subset_samples(wine.a, treatment == "Bare ground") # select bare ground samples
core.BG <- wine.core <- core(BG, detection = .0005, prevalence = .90) # a minimum abundance of 0.05 % and prevalent in 95% of the samples
core.BG = prune_taxa(taxa_sums(core.BG) > 0, core.BG) # only keep OTUs lager >0
core.BG <- row.names(tax_table(core.BG))

# Alternating cover 
AC <- subset_samples(wine.a, treatment == "Alternating cover") # select bare ground samples
core.AC <- wine.core <- core(AC, detection = .0005, prevalence = .90) # a minimum abundance of 0.05 % and prevalent in 95% of the samples
core.AC = prune_taxa(taxa_sums(core.AC) > 0, core.AC) # only keep OTUs lager >0
core.AC <- row.names(tax_table(core.AC))

# Complete cover 
CC <- subset_samples(wine.a, treatment == "Complete cover") # select bare ground samples
core.CC <- wine.core <- core(CC, detection = .0005, prevalence = .90) # a minimum abundance of 0.05 % and prevalent in 95% of the samples
core.CC = prune_taxa(taxa_sums(core.CC) > 0, core.CC) # only keep OTUs lager >0
core.CC <- row.names(tax_table(core.CC))

core=list(core.BG, core.AC, core.CC)

# Chart
venn.diagram(
  core,
  category.names = c("Bare ground" , "Alternating cover" , "Complete cover"),
  filename = 'Venn_corepertrtm.png',
  output=TRUE
)

####
venn.diagram(
  core,
  category.names = c("Bare ground" , "Alternating cover" , "Complete cover"),
  filename = 'Venn_coreptrtm_col.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 800 , 
  width = 800 , 
  resolution = 300,
  compression = "lzw",
  
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
  rotation = 1
)
