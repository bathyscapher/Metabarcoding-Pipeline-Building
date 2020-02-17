Remove the machine identifier & co. from the FASTQ filenames:
```bash
rename -n 's/MI.M03555_0\d{3}.001.FLD0\d{3}.SCI0\d{5}_//' *fastq.gz
rename -n 's/_17_16S//' *fastq.gz
```


