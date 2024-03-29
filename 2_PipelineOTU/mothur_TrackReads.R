################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### mothur pipeline OTU: Count read loss along pipeline
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2021
### korn@cumulonimbus.at
################################################################################
################################################################################
### To create countReadsMothur.csv, run in directory with mothur results:
### awk -f transposeList2Table.awk <(grep 'mothur > \|# of' mothur.*.logfile | grep '# of\|make.contigs\|screen.seqs\|align.seqs\|filter.seqs\|pre.cluster\|chimera.vsearch\|awk.\+single.accnos' | sed -r 's/mothur > /mothur: /' | sed -r 's/mothur/\nmothur/' | sed -r 's/\t/ /' | sed -r 's/#/Number/') > countReadsMothur.csv
################################################################################
################################################################################


library("ggplot2")
theme_set(theme_bw(base_size = 12))
library("reshape2")


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
reads <- read.table("countReads.csv", sep = "\t", header = TRUE)
reads <- reads[-1, ] # delete empty line
reads$mothur <- as.factor(reads$mothur)


levels(reads$mothur)[7] <- gsub("screen.seqs", "screen.seqs2", reads[4, 1])
reads$mothur <- gsub("\\(.+\\)", "", reads$mothur)


reads$mothur <- ordered(reads$mothur, levels = c("make.contigs", "screen.seqs",
                                                 "align.seqs", "screen.seqs2",
                                                 "filter.seqs", "pre.cluster",
                                                 "chimera.vsearch", "system"))


reads.m <- melt(reads, id.vars = "mothur", variable.name = "Counts",
                value.name = "Reads")


################################################################################
## Plot
ggplot(data = reads.m, aes(x = mothur, y = log(Reads), color = Counts)) +
  geom_line(aes(group = Counts), linetype = "dashed") +
  geom_point() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("")


################################################################################
################################################################################
