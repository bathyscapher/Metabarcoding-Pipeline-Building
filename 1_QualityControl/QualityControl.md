# Preprocessing & quality control
## Preprocess FASTQ
Remove the machine identifier & co. from the FASTQ filenames to keep only the actual sample names.
```bash
rename -n 's/MI.M03555_0\d{3}.001.FLD0\d{3}.SCI0\d{5}_//' *fastq.gz
rename -n 's/_17_16S//' *fastq.gz # for 16S
rename -n 's/_17_18S//' *fastq.gz # for 18S
```

## Filter & trim
Filter  and trim the reads with [this dada2 script](QualityFiltering.R).


