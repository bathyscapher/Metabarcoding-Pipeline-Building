# Customize SILVA db as reference alignment
For each primer pair, a customized SILVA SSU database 138.1 (`silva.full_v138_1.fasta`) was created following the instructions [here](https://mothur.org/blog/2021/SILVA-v138_1-reference-files/) and [here](https://mothur.org/blog/2016/Customization-for-your-region/).


Download the db and process it with ARB:
```
wget -N https://www.arb-silva.de/fileadmin/arb_web_db/release_138_1/ARB_files/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz
gunzip SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz
arb SILVA_138.1_SSURef_NR99_12_06_20_opt.arb
```

You may need to install [ARB](http://www.arb-home.de/) before:
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
Herein, the final output is `silva.v138.1_16S-V4.full` to classify the sequences.

The *E. coli* 16S rRNA gene ([NR_024570.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_024570.1/)) was trimmed to the amplified region by the respective primers (i.e. the sequence flanked by the forward and reverse primer). The primers need to be in an [`oligos` file](https://mothur.org/wiki/oligos_file/). For example, the Earth Microbiome Project 16S primers would look like in an `oligos` file as follows:
```
forward GTGYCAGCMGCCGCGGTAA
reverse ATTAGANACCCNNGTAGTCC
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



1. From the last command `summary.seqs` enter the coordinates as start and end position and run in `mothur`,
1. identify and 
1. grab the unique sequences without regard to their alignment,
1. rename fasta to `align` file and
1. extract the taxonomic information from the header.

```
screen.seqs(fasta=../silva.full_v138_1.fasta, start=13862, end=23444, maxambig=5)
pcr.seqs(start=13862, end=23444, keepdots=F)
degap.seqs()
unique.seqs()

system(grep ">" ../silva.full_v138_1.good.pcr.ng.unique.fasta | cut -f 1 | cut -c 2- > ../silva.full_v138_1.good.pcr.ng.unique.accnos)
get.seqs(fasta=../silva.full_v138_1.good.pcr.fasta, accnos=../silva.full_v138_1.good.pcr.ng.unique.accnos)
system(mv ../silva.full_v138_1.good.pcr.pick.fasta ../silva.v138_1_16S-V4.align)

## Refuses to write the file, thus, directly created in bash
grep '>' ../silva.v138_1_16S-V4.align | cut -f1,3 | cut -f2 -d'>' > silva.v138_1_16S-V4.full
```



## Build the SEED references
```
grep ">" silva.nr_v138.align | cut -f 1,2 | grep "\t100" | cut -f 1 | cut -c 2- > silva.seed_v138.accnos
mothur "#get.seqs(fasta=silva.nr_v138.align, taxonomy=silva.full_v138.tax, accnos=silva.seed_v138.accnos)"
mv silva.nr_v138.pick.align silva.seed_v138.align
mv ../silva.full_v138.pick.tax ../silva.seed_v138.tax

mothur "#get.seqs(taxonomy=silva.full_v138.tax, accnos=silva.full_v138.good.pcr.ng.unique.accnos)"
mv silva.full_v138.pick.tax silva.nr_v138.tax
```

