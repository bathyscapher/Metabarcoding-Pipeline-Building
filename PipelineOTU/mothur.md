# mothur pipeline
mothur runs directly in the shell (command line/terminal). Calling `mothur` opens the program and the commands below can be run from there.

Set working directory and specify number of available processors.
```
set.dir(tempdefault=.)
set.current(processors=6)
```

## Merge reads to contigs
Give a prefix for the output files and merge the forward and reverse reads to contigs.
```
make.file(inputdir=., type=gz, prefix=wine)
make.contigs(file=wine.files)
count.groups(group=wine.contigs.groups)
summary.seqs()
get.current()
```

## Quality trimming
Remove all contigs that are either too long or too short and exceed a certain amount of homopolymers. Length of EMP 16S amplicons = 300 to 350 bp.
```
screen.seqs(fasta=wine.trim.contigs.fasta, group=wine.contigs.groups, summary=wine.trim.contigs.summary, maxambig=0, minlength=250, maxlength=295, maxhomop=8)
count.groups(group=wine.contigs.good.groups)
summary.seqs(fasta=current)
get.current()
```

Count redundant sequences to speed up computation.
```
unique.seqs(fasta=wine.trim.contigs.good.fasta)
count.seqs(name=wine.trim.contigs.good.names, group=wine.contigs.good.groups)
summary.seqs(fasta=current, count=current)
```

## Align the contigs
Align the contigs to the (customized) reference data base (SILVA).
```
align.seqs(fasta=wine.trim.contigs.good.unique.fasta, reference=silva.v132_16S-V4.align, flip=f)
summary.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table)
get.current()
```

Remove all sequences from the alignment that lie outside ?? and update the count table and summary file.
```
screen.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table, summary=wine.trim.contigs.good.unique.summary, start=8, end=9581)
summary.seqs(fasta=current, count=current)
```

## Filter ...
```
filter.seqs(fasta=wine.trim.contigs.good.unique.good.align, vertical=T, trump=.)
unique.seqs(fasta=wine.trim.contigs.good.unique.good.filter.fasta, count=wine.trim.contigs.good.good.count_table)
summary.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.fasta, count=wine.trim.contigs.good.unique.good.filter.count_table)
system(rm -f wine.filter)
```

## Precluster sequences
Unify similar sequences (clusters sequences that are likely to be noisy!). A diff=2 refers to 4 nt for 2 seqs and thus, for 290 nt about 1.4 %. If diff = 3, 2.1 %, diff = 5, 3.4 %, diff
= 6, 4.1 %.
```
pre.cluster(fasta=wine.trim.contigs.good.unique.good.filter.unique.fasta, count=wine.trim.contigs.good.unique.good.filter.count_table, diffs=3, processors=6)
summary.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.count_table)
get.current()
```

## Detect chimeras
Detect chimeras with [vsearch](https://github.com/torognes/vsearch) and remove them.
```
chimera.vsearch(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
get.current()
remove.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
count.groups(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
summary.seqs(fasta=current, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
get.current()
```

## Remove singletons
Remove all reads that occur only once.
```
system(awk -e '{if ($2 < 50) { print $0 } }' wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)
remove.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)
#remove.seqs(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)
system(bash removeSeqsFromCountTable16S.sh)

count.groups(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
summary.seqs(fasta=current, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
get.current()
```

## Classify sequences
```
classify.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, reference=~/Desktop/SMP_unsynced/silva/16S-V4/silva.v132_16S-V4.align, taxonomy=~/Desktop/SMP_unsynced/silva/16S-V4/silva.v132_16S-V4.tax, cutoff=80)
get.current()
```

## Get phylotypes
```
phylotype(taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.taxonomy)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.list, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.org.list)
get.current()
```

## Extract OTUs at different taxa level
Remove the multiple in-line headers from file. FNR: the record number (typically the line number).
```
system(awk -e '{if (($1 != "label") || (FNR == 1)) { print $0 } }' wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.org.list > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.list)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=1)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.1.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=2)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.2.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=3)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.3.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=4)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.4.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=5)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.5.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=6)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.tx.6.shared)

get.current()
```

## Classify OTUs
Find consensus taxonomy.
```
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.taxonomy)
get.current()
```

## Compute distance matrix and cluster OTUs
Note: avoid `cluster.split`, results are not accurate.
```
dist.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.04)
get.current()

cluster(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03, method=opti)
get.current()
```

## Get OTU table
Set the threshold (by convention a 3 % barcoding gap) and create the OTU-table.
```
make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#summary.single(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)
```

## Count OTUs per sample
**awk:** -v: var=val, -e: use program-text, OFS: output field separator, NF: input field number, NR: total number of input records so far.
```
system(awk -v OFS='\t' -e '{notus=0; for (i=4; i<=NF; i++) { if ($i > 0) notus++; }; if (NR > 1) print $2 OFS notus; }' wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.sample.summary)
get.current()
```

## Get representative sequences for OTUs
= cluster centroids?
```
get.oturep(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
get.oturep(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
```

## Classify clustered OTUs
```
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_16S_V4.wang.taxonomy, label=0.03)
get.current()
```

## Rarefaction curves
Rarefaction curves show if the sampling (here, sequencing) depth is sufficient.
```
rarefaction.single(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, calc=sobs, freq=100)
```
