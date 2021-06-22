################################################################################
################################################################################
################################################################################
################################################################################
### Source: http://blog.mothur.org/2018/01/10/SILVA-v132-reference-files/
### Map the taxa from SILVA to the six Linnean levels (kingdom, phylum, class,
### order, family, and genus).
### Author: Eric Collins, University of Alaska Fairbanks
### Modified: korn@cumulonimbus.at 2019, 2020, 2021


rm(list = ls())
setwd("~/docker/silva/")


map.in <- read.table("tax_slv_ssu_138.1.txt", header = FALSE, sep = "\t",
           stringsAsFactors = FALSE)
map.in <- map.in[, c(1, 3)]
colnames(map.in) <- c("Taxa", "Rank")


ranks <- c("root", "domain", "major_clade", "superkingdom", "kingdom",
           "subkingdom", "infrakingdom", "superphylum", "phylum", "subphylum",
           "infraphylum", "superclass", "class", "subclass", "infraclass",
           "superorder", "order", "suborder", "superfamily", "family",
           "subfamily", "genus")

rankAbb <- c("ro", "do", "mc", "pk", "ki", "bk", "ik", "pp", "ph", "bp", "ip",
             "pc", "cl", "bc", "ic", "po", "or", "bo", "pf", "fa", "bf", "ge")

tax.mat <- matrix(data = "",
                  nrow = nrow(map.in),
                  ncol = length(ranks))
tax.mat[, 1] <- "root"
colnames(tax.mat) <- ranks


outlevels <- c("domain", "phylum", "class", "order", "family", "genus")


for(i in 1:nrow(map.in)) {
  taxname <- unlist(strsplit(as.character(map.in[i, 1]), split = ';'))
  #print(taxname);
  while (length(taxname) > 0) {
    # regex to look for exact match
    tax.exp <- paste(paste(taxname, collapse = ";"), ";", sep = "")
    tax.match <- match(tax.exp, map.in$Taxa)
    tax.mat[i, map.in[tax.match, 2]] <- tail(taxname, 1)
    taxname <- head(taxname, -1)
    }
}



## Fill in the gaps by using the closest higher taxonomic level appended with
## an abbreviation for the current taxonomic level (skip, if not needed)
for(i in 1:nrow(tax.mat)) {
  for(j in 1:ncol(tax.mat)) {
    if(tax.mat[i, j] < 0) {
      tax.mat[i, j] <- paste(tmptax, rankAbb[j], sep = "_")
      }
    else {
      tmptax <- tax.mat[i, j]
      }
    }
  ## Map the new name to the input taxonomic levels
  map.in[i, "taxout"] <- paste(paste(tax.mat[i, outlevels], collapse = ";"),
                               ";", sep = "")
  }


## Replace spaces with underscores
map.in$taxout <- gsub(" ", "_", map.in$taxout)


################################################################################
## Import the taxonomic levels from SILVA and remap them with the new levels
## 16S
tax.in <- read.table("prok/silva.v138.1_16S-V4.full", header = FALSE,
                     stringsAsFactors = FALSE, sep = "\t")

## 18S
tax.in <- read.table("euk/silva.v138.1_18S-V4.full", header = FALSE,
                     stringsAsFactors = FALSE, sep = "\t")

colnames(tax.in) <- c("taxaID", "Taxa")


## Correct the "...;Polaribacter;Polaribacter;" problem
# tax.in$Taxa <- gsub("Polaribacter;Polaribacter;", "Polaribacter;",
#                         tax.in$Taxa)
# tax.in$Taxa <- gsub("Polaribacter;Polaribacter 3;", "Polaribacter 3;",
#                         tax.in$Taxa)
tax.in$Taxa <- gsub(";[[:space:]]+$", ";", tax.in$Taxa)


## Replace double semicolon
tax.in$Taxa <- gsub("Bacteria;GBS-1;;", "Bacteria;GBS-1;", tax.in$Taxa)

## Insert missing final ';'
tax.in$Taxa <- gsub("D64120.DX3Polym", "D64120.DX3Polym;", tax.in$Taxa)

## Remove white space in 'Incertae sedis'
tax.in$Taxa <- gsub("Incertae sedis;", "Incertae_Sedis;", tax.in$Taxa)


tax.in$id <- 1:nrow(tax.in)

tax.write <- merge(tax.in, map.in, all.x = TRUE, sort = FALSE)
tax.write <- tax.write[order(tax.write$id), ]

tax.write[is.na(map.in$taxout),]


# Check if everything has 6 taxonomic levels (kingdom to genus)
getDepth <- function(taxonString){
  initial <- nchar(taxonString)
  removed <- nchar(gsub(";", "", taxonString))
  return(initial - removed)
}



depth <- getDepth(tax.write$taxout)
summary(depth) # should all be 6 and there should be no NAs


tax.write[is.na(tax.write$taxout), ]



bacteria <- grepl("Bacteria;", tax.write$taxout)
archaea <- grepl("Archaea;", tax.write$taxout)
eukarya <- grepl("Eukaryota;", tax.write$taxout)


tax.write[depth > 6 & bacteria,] # good to go
tax.write[depth > 6 & archaea,] # good to go
tax.write[depth > 6 & eukarya,] # good to go


## 16S V4
write.table(tax.write[, c("taxaID", "taxout")],
            file = "prok/silva.v138.1_16S-V4.tax",
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

## 18S V4
write.table(tax.write[, c("taxaID", "taxout")],
            file = "euk/silva.v138.1_18S-V4.tax",
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


################################################################################
################################################################################
