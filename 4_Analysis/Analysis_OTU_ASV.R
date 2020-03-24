################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### mothur pipeline OTU: rarefaction curves
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2020
### korn@cumulonimbus.at
################################################################################
################################################################################

rm(list = ls())


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))


setwd("/scratch/PromESSinG/")

################################################################################
### OTU

list.files(pattern = "3.cons.taxo", recursive = TRUE)
# show_mothur_cutoffs("pilotFilt18S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list")

wine <- import_mothur(mothur_list_file = NULL,
                      mothur_group_file = NULL,
                      mothur_tree_file = NULL,
                      cutoff = NULL,
                      mothur_shared_file = "prok/filtered/wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared",
                      mothur_constaxonomy_file = "prok/filtered/wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy",
                      parseFunction = parse_taxonomy_default)


################################################################################
### ASV
seqtab.nochim <- readRDS("prok/seq.tab.nochim.rds")
tax.id <- readRDS("prok/tax.id.rds")


wine <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                 tax_table(tax.id))


dna <- Biostrings::DNAStringSet(taxa_names(wine))
names(dna) <- taxa_names(wine)
wine <- merge_phyloseq(wine, dna)
taxa_names(wine) <- paste0("ASV", seq(ntaxa(wine)))


################################################################################
### ...




################################################################################
################################################################################

