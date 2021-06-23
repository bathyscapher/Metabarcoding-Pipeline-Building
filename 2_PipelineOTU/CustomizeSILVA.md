# Customize SILVA db as reference alignment
For each primer pair, a customized SILVA SSU database 138.1 (`silva.full_v138.1.fasta`) was created following the instructions [here](https://mothur.org/blog/2021/SILVA-v138_1-reference-files/) and [here](https://mothur.org/blog/2016/Customization-for-your-region/).


Download the db and process it with [ARB](http://www.arb-home.de/):
```
wget -N https://www.arb-silva.de/fileadmin/arb_web_db/release_138.1/ARB_files/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz
gunzip SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz
arb SILVA_138.1_SSURef_NR99_12_06_20_opt.arb
```

You may need to install ARB before:
```
sudo apt update
sudo apt install arb
```

In ARB,

1. Click the search button
1. Set the first search field to `ARB_color` and set it to `1`. Click on the equal sign until it indicates inequality (this removes low quality reads and chimeras).
1. Click `Search`.
1. Click the `Mark Listed Unmark Rest` button.
1. Close the `Search and Query` box.
1. Click on `File > export > export to external format`. Set the `Export` option to `marked`, `Filter` to `none` and `Compress` to `no`. In the field for `Choose an output file name` enter your path and name the output file (e.g. `silva.full_v138.1.fasta`).
1. Create a custom formatting file that includes the sequences accession number and it's taxonomy across the top line: create `fasta_mothur.eft` in `/usr/lib/arb/lib/export` with the following content (depending on your OS, the path might deviate, on Windows it might be something like `$ARBHOME/lib/export/`):
    ```
    SUFFIX          fasta
    BEGIN
    >*(acc).*(name)\t*(align_ident_slv)\t*(tax_slv);
    *(|export_sequence)
    ```
1. In `Select a format` select `fasta_mothur.eft`.
1. Press `Go`.
1. Quit ARB.


## 16S
Herein, the final output is `silva.v138.1_16S-V4.align` to classify the sequences.

The *E. coli* 16S rRNA gene ([NR_024570.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_024570.1/)) was trimmed to the amplified region by the respective primers (i.e. the sequence flanked by the forward and reverse primer). The primers need to be in an [`oligos` file](https://mothur.org/wiki/oligos_file/). For example, the [Earth Microbiome Project 16S primers](https://earthmicrobiome.org/protocols-and-standards/16s/) would look in an `oligos` file as follows:
```
forward GTGYCAGCMGCCGCGGTAA
reverse ATTAGANACCCNNGTAGTCC
```
or:
```
forward GTGYCAGCMGCCGCGGTAA ATTAGANACCCNNGTAGTCC
```

1. Crop the *E. coli* fasta with the 16S-EMP primers,
1. Align it to [`silva.seed_v138.align`](https://mothur.org/wiki/silva_reference_files/) and extract relevant columns.

In `mothur`, run:
```
set.current(processors=6)
set.dir(tempdefault=.)
pcr.seqs(fasta=e-coli.fasta, oligos=EMP_16S.oligos)
system(mv e-coli.pcr.fasta e.coli.V4.fasta)
align.seqs(fasta=e.coli.V4.fasta, reference=../silva.seed_v138_1.align)
summary.seqs(fasta=e.coli.V4.align)
```


Run in `mothur`:
1. Entering the coordinates from the last `summary.seqs` command as start and end position removes sequences that are not full length (additionally, filter for sequences with more than 5 ambiguous base calls).
1. `pcr.seqs` converts all bases that occur before the given start- and after the given stop position to `.`, so that the sequences only span the region between the given primer binding sites.
1. Unalign the sequences (`degap.seqs`) and
1. identify the unique sequences (`unique.seqs`).
1. Convert the resulting `fasta` file into an `accnos` file
1. what allows to pull out the unique sequences from the aligned file (`get.seqs`),
1. rename the `fasta` to `align` file and
1. extract the taxonomic information from the header (for some reason, `mothur 1.45.3` refuses to write the file, thus, directly create it in bash (should work in the next `mothur` version)):
```
# in mothur
screen.seqs(fasta=../silva.full_v138.1.fasta, start=13862, end=23444, maxambig=5)
pcr.seqs(start=13862, end=23444, keepdots=F)
degap.seqs()
unique.seqs()

system(grep ">" ../silva.full_v138.1.good.pcr.ng.unique.fasta | cut -f 1 | cut -c 2- > ../silva.full_v138.1.good.pcr.ng.unique.accnos)
get.seqs(fasta=../silva.full_v138.1.good.pcr.fasta, accnos=../silva.full_v138.1.good.pcr.ng.unique.accnos)
system(command='mv ../silva.full_v138.1.good.pcr.pick.fasta ../silva.v138.1_16S-V4.align')

# in bash
grep '>' ../silva.v138.1_16S-V4.align | cut -f1,3 | cut -f2 -d'>' > silva.v138.1_16S-V4.full
```


## 18S
Herein, the final output is `silva.v138.1_18S-V4.align` to classify the sequences.

The very same as for the 16S above, but with the *S. cerevisae* 18S rRNA gene ([NR_132222.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_132222.1?report=fasta)) and the [Earth Microbiome Project 18S primers](https://earthmicrobiome.org/protocols-and-standards/18s/):
```
forward GTACACACCGCCCGTC
reverse TGATCCTTCYGCAGGTTCACCTAC
```
or:
```
forward GTACACACCGCCCGTC TGATCCTTCYGCAGGTTCACCTAC
```


In `mothur` run:
```
set.current(processors=6)
set.dir(tempdefault=.)
pcr.seqs(fasta=s-cerevisiae.fasta, oligos=EMP_18S.oligos)
system(command='mv s-cerevisiae.pcr.fasta s.cerevisiae.V4.fasta')
align.seqs(fasta=s.cerevisiae.V4.fasta, reference=../silva.seed_v138_1.align)
summary.seqs(fasta=s.cerevisiae.V4.align)


screen.seqs(fasta=../silva.full_v138.1.fasta, start=42554, end=43116, maxambig=5)
pcr.seqs(start=42554, end=43116, keepdots=F)
degap.seqs()
unique.seqs()

system(command='grep ">" ../silva.full_v138.1.good.pcr.ng.unique.fasta | cut -f 1 | cut -c 2- > ../silva.full_v138.1.good.pcr.ng.unique.accnos')
get.seqs(fasta=../silva.full_v138.1.good.pcr.fasta, accnos=../silva.full_v138.1.good.pcr.ng.unique.accnos)
system(command='mv ../silva.full_v138.1.good.pcr.pick.fasta ../silva.v138.1_18S-V4.align')
summary.seqs(fasta=../silva.v138.1_18S-V4.align)
```

For some reason, `mothur 1.45.3` refuses to write the file, thus, directly create it in bash (should work in the next `mothur` version):
```
grep '>' ../silva.v138.1_18S-V4.align | cut -f1,3 | cut -f2 -d'>' > ../silva.v138.1_18S-V4.full
```


## Formatting the taxonomy files
Herein, the final output is `silva.v138.1_16S-V4.tax` and `silva.v138.1_18S-V4.tax`, respectively, to taxonomically classify the sequences.


First, download the SILVA taxa mapping file:
```
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz
gunzip tax_slv_ssu_138.1.txt.gz
```

To format the taxonomy file so that it is compatible with `mothur`, run in `R`:
```
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
```


