Remove the machine identifier & co. from the FASTQ filenames:
```bash
rename -n 's/MI.M03555_0\d{3}.001.FLD0\d{3}.SCI0\d{5}_//' *fastq.gz
rename -n 's/_17_16S//' *fastq.gz
```


Filter  and trim the reads. Place filtered FASTQ in `filtered` subdirectory. Exemplarily, compare read quality profiles raw and filtered.
```R
rF.f <- file.path("filtered", paste0(sample.names, "_16S_R1_filt.fastq.gz"))
rR.f <- file.path("filtered", paste0(sample.names, "_16S_R2_filt.fastq.gz"))

names(rF.f) <- sample.names
names(rR.f) <- sample.names

out <- filterAndTrim(rF, rF.f, rR, rR.f,
                     maxN = 0, maxEE = c(2,2), truncQ = 10, rm.phix = TRUE,
                     compress = TRUE, multithread = ncore)

plotQualityProfile(rF[1:3])
plotQualityProfile(rF.f[1:3])
plotQualityProfile(rR[1:3])
plotQualityProfile(rR.f[1:3])
```



