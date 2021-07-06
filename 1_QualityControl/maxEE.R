################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### Calculate and plot maxEE
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2021
### korn@cumulonimbus.at
################################################################################
################################################################################


library("ShortRead")


rm(list = ls())


## Choose pro- or eukaryotes
primer <- "16S"
# primer <- "18S"


if (primer == "16S") {
  setwd("/home/rstudio/prok/")
} else {setwd("/home/rstudio/euk/")}
getwd()


################################################################################
################################################################################

calcEvalue <- function(qvalue) {
  evalue <- 0
  for (i in 1:length(qvalue)) {
    if (!is.na(qvalue[i])) {
      evalue <- evalue + 10^(-qvalue[i]/10)
    }
  }
  return(evalue)
}
calcEvalue(c(30,31,34,23,10, NA))
calcEvalue(c(3,5,10,20,10, NA))


getEvalues <- function(fname)
{
  my.fastq <- readFastq(fname)
  my.fastq.q <- quality(my.fastq)
  my.fastq.scores <- as(my.fastq.q, "matrix")
  evalues <- apply(my.fastq.scores, 1, calcEvalue)
  return(evalues)
}

filtered.e <- getEvalues("filtered/101_16S_R1_filt.fastq.gz"); hist(filtered.e)
unfiltered.e <- getEvalues("filtN/MI.M03555_0244.001.FLD0157.SCI019503_101_17_16S_R1_sub.fastq.gz")
par(mfrow = c(1, 2))
hist(unfiltered.e)


my.fastq.scores <- as(my.fastq.q, "matrix")
head(my.fastq.scores)

calcEvalue(my.fastq.scores[1,])
my.fastq.scores[1,]


my.fastq.e <- apply(my.fastq.scores, 1, calcEvalue)
head(my.fastq.e)

hist(my.fastq.e)

################################################################################
################################################################################