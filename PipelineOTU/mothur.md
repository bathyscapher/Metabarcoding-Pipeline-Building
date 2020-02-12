# mothur pipeline
mothur runs directly in the shell (command line/terminal). Calling `mothur` opens the program and the commands below can be run from there.

Set working directory and specify number of available processors.
```bash
set.dir(tempdefault=.)
set.current(processors=6)
```

## Merge reads to contigs
Create a metafile listing all the fastq in the directory with a given prefix (`wine.files`). Then, merge the forward and reverse reads to contigs with the [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) and show summary statistics.
```bash
make.file(inputdir=., type=gz, prefix=wine)
make.contigs(file=wine.files, align=needleman)
count.groups(group=wine.contigs.groups)
summary.seqs()
get.current()
```

## Quality trimming
Remove all contigs that are either too long or too short and exceed a certain amount of homopolymers. Values taken from previous `summary.seqs()`. ???Amplicon length of the [EMP 16S primers](http://www.earthmicrobiome.org/protocols-and-standards/16s/) is about 390 bp???
```bash
screen.seqs(fasta=wine.trim.contigs.fasta, group=wine.contigs.groups, summary=wine.trim.contigs.summary, maxambig=0, minlength=250, maxlength=295, maxhomop=8)
count.groups(group=wine.contigs.good.groups)
summary.seqs(fasta=current)
get.current()
```

Count and exclude redundant sequences to speed up computation.
```bash
unique.seqs(fasta=wine.trim.contigs.good.fasta)
count.seqs(name=wine.trim.contigs.good.names, group=wine.contigs.good.groups)
summary.seqs(fasta=current, count=current)
```

## Align the contigs
Align the contigs to the (customized) SILVA reference data base.
```bash
align.seqs(fasta=wine.trim.contigs.good.unique.fasta, reference=silva.v132.align, flip=f)
summary.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table)
get.current()
```

Remove all sequences from the alignment that are outliers??.
```bash
screen.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table, summary=wine.trim.contigs.good.unique.summary, start=8, end=9581)
summary.seqs(fasta=current, count=current)
```

## Remove empty columns
Remove columns that contain only gaps and grap unique sequences only.
```bash
filter.seqs(fasta=wine.trim.contigs.good.unique.good.align, vertical=T, trump=.)
unique.seqs(fasta=wine.trim.contigs.good.unique.good.filter.fasta, count=wine.trim.contigs.good.good.count_table)
summary.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.fasta, count=wine.trim.contigs.good.unique.good.filter.count_table)
system(rm -f wine.filter)
```

## Precluster sequences
Unify similar sequences, i.e. cluster sequences that probably are noisy.

**diff** | **2** | **3** | **5** | **6**
:--- | ---: | ---: | ---: | ---:
**nt for 2 seqs** | 4 | 6 | 10 | 12 
**% for 290 nt** | 1.4 | 2.1 | 3.4 | 4.1 

```bash
pre.cluster(fasta=wine.trim.contigs.good.unique.good.filter.unique.fasta, count=wine.trim.contigs.good.unique.good.filter.count_table, diffs=3, processors=6)
summary.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.count_table)
get.current()
```

## Detect chimeras
Detect chimeras with [vsearch](https://github.com/torognes/vsearch) and remove them.
```bash
chimera.vsearch(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
get.current()
remove.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
count.groups(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
summary.seqs(fasta=current, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
get.current()
```

## Remove singletons
Remove all reads that occur only once and update statistics.
```bash
system(awk -e '{if ($2 < 50) { print $0 } }' wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)

remove.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)
remove.seqs(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)
#system(bash removeSeqsFromCountTable16S.sh)

count.groups(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
summary.seqs(fasta=current, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
get.current()
```

## Classify sequences
Assign sequences to the reference database and taxonomy with the [Ribosomal Database Project (RDP) classifier](https://aem.asm.org/content/73/16/5261).
```bash
classify.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, reference=silva.v132.align, taxonomy=silva.v132.tax, method=wang, cutoff=80)
get.current()
```

## Get phylotypes
Assign sequences to OTUs based on their taxonomy and outputs a .list, .rabund and .sabund file.
```bash
phylotype(taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.taxonomy)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.list, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.org.list)
get.current()
```

## Extract OTUs at different taxa level
First, remove the multiple in-line headers from the file. Then, create `.shared` files for all six taxonomic levels (1 = , 2 = , 3 = , 4 = , 5 = , 6 = ).
```bash
system(awk -e '{if (($1 != "label") || (FNR == 1)) { print $0 } }' wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.org.list > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.list)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=1)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.1.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=2)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.2.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=3)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.3.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=4)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.4.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=5)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.5.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=6)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.tx.6.shared)

get.current()
```

## Classify OTUs
Find consensus taxonomy for an OTU.
```bash
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.taxonomy)
get.current()
```

## Compute distance matrix and cluster OTUs
Calculate uncorrected pairwise distances between the aligned DNA sequences; distances > 0.04 are ignored to safe resources (note: avoid `cluster.split`, results are not accurate). Then, assign this sequences to OTUs with clustering using the ??`opti`?? method.
```bash
dist.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.04)
get.current()

cluster(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03, method=opti)
get.current()
```

## Get OTU table
Create a .shared and a .rabund file for each group. Set the threshold (by convention a 3 % barcoding gap) and create the OTU-table.
```bash
make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#summary.single(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)
```

## Count OTUs per sample
Use [awk](https://en.wikipedia.org/wiki/AWK) to count &alpha;-diversity sample-wise (awk -v: var=val, -e: use program-text, OFS: output field separator, NF: input field number, NR: total number of input records so far).
```bash
system(awk -v OFS='\t' -e '{notus=0; for (i=4; i<=NF; i++) { if ($i > 0) notus++; }; if (NR > 1) print $2 OFS notus; }' wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.sample.summary)
get.current()
```

## Get representative sequences for OTUs
= cluster centroids?

Returns a FASTA file where the headers additionally carry the OTU number and the total abundance.
```bash
#get.oturep(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

get.oturep(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
```

## Classify clustered OTUs
Get consenus taxonomy for OTUs.
```bash
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.taxonomy, label=0.03)
get.current()
```

## Rarefaction curves
Calculate sample-wise rarefaction curves to control for sufficient sampling (here, sequencing) depth.
```bash
rarefaction.single(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, calc=sobs, freq=100)
```
