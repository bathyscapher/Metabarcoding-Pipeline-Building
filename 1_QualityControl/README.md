# Preprocessing & quality control
## Preprocess FASTQ
Remove the machine identifier & co. in bash from the FASTQ filenames to keep only the actual sample names. In the directory with the files run:
```bash
rename -n 's/MI.M03555_0\d{3}.001.FLD0\d{3}.SCI0\d{5}_//' *fastq.gz
rename -n 's/_17_1\d{1}S//' *fastq.gz
```

## Filter & trim the reads
Filter  and trim the reads with [this R script](QualityFiltering.R).

![16S reads unfiltered vs. filtered](/Graphs/Preprocessing_QualityFiltering16S.png)
