################################################################################
################################################################################
### Metabarcoding Pipeline Building: CUSO Workshop
### mothur pipeline OTU: rarefaction curves
### Gerhard Thallinger, Rachel Korn & Magdalena Steiner 2020
### korn@cumulonimbus.at
################################################################################
################################################################################

rm(list = ls())

library("ggplot2")
theme_set(theme_bw(base_size = 12)
library("reshape2")

setwd("~")

################################################################################
# list.files(path = ".", pattern = "rare", recursive = TRUE)

### 16S
rarefied16S <- read.table("16S/wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.rarefaction"
                          ,
                          sep = "\t", header = TRUE, check.names = FALSE)

numsampled <- rarefied16S$numsampled

rarefied16S.sample <- cbind(numsampled,
                            rarefied16S[ , grepl("0.03", names(rarefied16S))])
rarefied16S.lci <- cbind(numsampled,
                         rarefied16S[ , grepl("lci", names(rarefied16S))])
rarefied16S.hci <- cbind(numsampled,
                         rarefied16S[ , grepl("hci", names(rarefied16S))])


rarefied16S.sample.m <- melt(rarefied16S.sample, id = "numsampled",
                             variable.name = "Sample",
                             value.name = "OTU")
rarefied16S.lci.m <- melt(rarefied16S.lci, id = "numsampled",
                          variable.name = "Sample",
                          value.name = "lci")
rarefied16S.hci.m <- melt(rarefied16S.hci, id = "numsampled",
                          variable.name = "Sample",
                          value.name = "hci")


rarefied16S.sample.m$Sample <- as.factor(gsub("0.03-", "",
                                              rarefied16S.sample.m$Sample))
rarefied16S.lci.m$Sample <- as.factor(gsub("lci-", "", rarefied16S.lci.m$Sample))
rarefied16S.hci.m$Sample <- as.factor(gsub("hci-", "", rarefied16S.hci.m$Sample))


rarefied16S.merged <- merge(rarefied16S.sample.m, rarefied16S.lci.m,
                            by = c("numsampled", "Sample"))
rarefied16S.merged <- merge(rarefied16S.merged, rarefied16S.hci.m,
                            by = c("numsampled", "Sample"))

rm(rarefied16S, rarefied16S.hci, rarefied16S.hci.m, rarefied16S.lci,
   rarefied16S.lci.m, rarefied16S.sample, rarefied16S.sample.m, numsampled)


rarefied16S.merged$Site <- as.factor(gsub(".{5}$", "",
                                          rarefied16S.merged$Sample))
# rarefied16S.merged$Site <- ordered(rarefied16S.merged$Site,
#                                    levels = c("CB", "LT", "LE", "LV", "LM",
#                                               "MC", "NC"))
rarefied16S.merged$Succession <- as.factor(gsub("^.{6}", "",
                                                rarefied16S.merged$Sample))
rarefied16S.merged$Primer <- as.factor("16S")


################################################################################
### 18S
rarefied18S <- read.table("18S/wine.filter.cat.unique.good.filter.unique.precluster.pick.pick.opti_mcc.groups.rarefaction",
                          sep = "\t", header = TRUE, check.names = FALSE)

numsampled <- rarefied18S$numsampled


rarefied18S.sample <- cbind(numsampled,
                            rarefied18S[ , grepl("0.03", names(rarefied18S))])
rarefied18S.lci <- cbind(numsampled,
                         rarefied18S[ , grepl("lci", names(rarefied18S))])
rarefied18S.hci <- cbind(numsampled,
                         rarefied18S[ , grepl("hci", names(rarefied18S))])



rarefied18S.sample.m <- melt(rarefied18S.sample, id = "numsampled",
                             variable.name = "Sample",
                             value.name = "OTUs")
rarefied18S.lci.m <- melt(rarefied18S.lci, id = "numsampled",
                          variable.name = "Sample",
                          value.name = "lci")
rarefied18S.hci.m <- melt(rarefied18S.hci, id = "numsampled",
                          variable.name = "Sample",
                          value.name = "hci")


rarefied18S.sample.m$Sample <- as.factor(gsub("0.03-", "",
                                              rarefied18S.sample.m$Sample))
rarefied18S.lci.m$Sample <- as.factor(gsub("lci-", "", rarefied18S.lci.m$Sample))
rarefied18S.hci.m$Sample <- as.factor(gsub("hci-", "", rarefied18S.hci.m$Sample))


rarefied18S.merged <- merge(rarefied18S.sample.m, rarefied18S.lci.m,
                            by = c("numsampled", "Sample"))
rarefied18S.merged <- merge(rarefied18S.merged, rarefied18S.hci.m,
                            by = c("numsampled", "Sample"))

rm(rarefied18S, rarefied18S.hci, rarefied18S.hci.m, rarefied18S.lci,
   rarefied18S.lci.m, rarefied18S.sample, rarefied18S.sample.m, numsampled)


rarefied18S.merged$Site <- as.factor(gsub(".{5}$", "",
                                          rarefied18S.merged$Sample))
# rarefied18S.merged$Site <- ordered(rarefied18S.merged$Site,
#                                    levels = c("CB", "LT", "LE", "LV", "LM",
#                                               "MC", "NC"))
rarefied18S.merged$Succession <- as.factor(gsub("^.{6}", "",
                                                rarefied18S.merged$Sample))

rarefied18S.merged$Primer <- as.factor("18S")

################################################################################
## Merge 16S and 18S

rarefied <- rbind(rarefied16S.merged, rarefied18S.merged)
head(rarefied)

levels(rarefied$Succession) <- c("Early succession", "Late succession", "Moss",
                                 "Negative control Pilot", "Negative control SMP",
                                 "Mock Prokaryote Pilot", "Mock Prokaryote SMP",
                                 "reseq Pilot",
                                 "Mock Eukaryote Pilot", "Mock Eukaryote SMP")

################################################################################
## Plot rarefaction curves

ggplot(data = rarefied16S.merged) +
  geom_line(aes(x = numsampled, y = lci, group = Sample), color = "gray85",
            size = 0.3) +
  geom_line(aes(x = numsampled, y = hci, group = Sample), color = "gray85",
            size = 0.3) +
  geom_line(aes(x = numsampled, y = OTUs, color = Succession, group = Sample),
            size = 0.5) +
  facet_grid(Primer ~ Site) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Reads") +
  ylab("OTUs")

################################################################################
################################################################################
