# mothur pipeline
## 16S & 18S
The pipeline is basically the same for 16S and 18S. Exceptions are the SILVA reference database as it was customized to each primer pair and some consequential steps - such cases are highlighted. 

Open `mothur` directly in the shell (command line/terminal). Calling `mothur` opens the program and the commands below can be run from there.

Set working directory and specify number of available processors.
```bash
set.dir(tempdefault=.)
set.current(processors=6)
```
*Note*: run `nproc` to see number of available processors.

## Merge reads to contigs
Create a metafile listing all the fastq in the directory with a given prefix (`wine.files`). Then, merge the forward and reverse reads to contigs with the [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).
```bash
make.file(inputdir=., type=gz, prefix=wine)
make.contigs(file=wine.files, align=needleman)

count.groups(group=wine.contigs.groups)
summary.seqs()
get.current()
```

## Quality trimming
Remove all contigs that are either too long or too short and exceed a certain amount of homopolymers. Values taken from previous `summary.seqs()`.

The expected amplicon length of the [EMP 16S primers](http://www.earthmicrobiome.org/protocols-and-standards/16s/) is about 300 to 350 bp, the [EMP 18S primers](http://www.earthmicrobiome.org/protocols-and-standards/18s/) about 260 +/- 50 bp.
```bash
### 16S
screen.seqs(fasta=wine.trim.contigs.fasta, group=wine.contigs.groups, summary=wine.trim.contigs.summary, maxambig=0, minlength=290, maxlength=295, maxhomop=8)

### 18S
screen.seqs(fasta=wine.trim.contigs.fasta, group=wine.contigs.groups, summary=wine.trim.contigs.summary, maxambig=0, minlength=300, maxlength=340, maxhomop=15)

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

## Align the contigs to SILVA
* How can the range in coordinates stretch so far over SILVA?

Align the contigs to the (customized) [SILVA](https://www.arb-silva.de/) reference data base.

*Note:* instructions to taylor the SILVA database to specific primers are [here](http://blog.mothur.org/2018/01/10/SILVA-v132-reference-files/) and [here](http://blog.mothur.org/2016/07/07/Customization-for-your-region/).

```bash
### 16S
align.seqs(fasta=wine.trim.contigs.good.unique.fasta, reference=silva.v132_EMP16S.align, flip=f)

### 18S
align.seqs(fasta=wine.trim.contigs.good.unique.fasta, reference=silva.v132_EMP18S.align, flip=f)

summary.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table)
get.current()
```

Screen once more, to ensure that all sequences overlap the same region.
```bash
### 16S
screen.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table, summary=wine.trim.contigs.good.unique.summary, start=8, end=9581)

### 18S
screen.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table, summary=wine.trim.contigs.good.unique.summary, start=4, end=562)

summary.seqs(fasta=current, count=current)
```

Remove columns that contain only gaps (they stem from the alignment) and grab unique sequences only.
```bash
filter.seqs(fasta=wine.trim.contigs.good.unique.good.align, vertical=T, trump=.)
unique.seqs(fasta=wine.trim.contigs.good.unique.good.filter.fasta, count=wine.trim.contigs.good.good.count_table)

summary.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.fasta, count=wine.trim.contigs.good.unique.good.filter.count_table)
system(rm -f wine.filter)
```

## Precluster sequences
Unify similar sequences, i.e. cluster sequences that probably are noisy aka "denoising" (according to the MiSeq SOP, allow 1 difference for every 100 bp of sequence).

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
remove.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

count.groups(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
summary.seqs(fasta=current, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
get.current()
```

## Remove singletons
Generate a file containing all singletons with [AWK](https://en.wikipedia.org/wiki/AWK) and remove them from the current fasta and count_table.
```bash
system(awk -e '{if ($2 < 2) {print $1}}' wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)

remove.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)
remove.seqs(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)

count.groups(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
summary.seqs(fasta=current, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
get.current()
```

## Classify sequences
* Only to remove lineages?

Assign sequences to the reference database and taxonomy with the [Ribosomal Database Project (RDP) classifier](https://aem.asm.org/content/73/16/5261).
```bash
### 16S
classify.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, reference=silva.v132_EMP16S.align, taxonomy=silva.v132_EMP16S.tax, method=wang, cutoff=80)

### 18S
classify.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, reference=silva.v132_EMP18S.align, taxonomy=silva.v132_EMP18S.tax, cutoff=80)

get.current()
```

## Get phylotypes
Assign sequences to OTUs based on their taxonomy and create a .list, .rabund and .sabund file.
```bash
### 16S
phylotype(taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.taxonomy)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.list, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.org.list)

### 16S
phylotype(taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.taxonomy)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.list, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.org.list)

get.current()
```

## Extract OTUs at different taxa level
* We could skip this, much simpler in phyloseq/R
First, remove the multiple in-line headers from the file with [AWK](https://en.wikipedia.org/wiki/AWK). Then, create .shared files for all six taxonomic levels (1 = domain, 2 = phylum, 3 = order, 4 = class?, 5 = family, 6 = genus).
```bash
### 16S
system(awk -e '{if (($1 != "label") || (FNR == 1)) {print $0}}' wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.org.list > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.list)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=1)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.1.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=2)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.2.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=3)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.3.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=4)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.4.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=5)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.5.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=6)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.tx.6.shared)

### 18S
system(awk -e '{if (($1 != "label") || (FNR == 1)) {print $0}}' wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.org.list > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.list)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=1)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.1.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=2)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.2.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=3)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.3.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=4)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.4.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=5)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.5.shared)

make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=6)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared)
rename.file(input=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.shared, new=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.tx.6.shared)

get.current()
```

## Classify OTUs
* Not only after creating the distance matrix?

Find consensus taxonomy for an OTU.
```bash
### 16S
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.taxonomy)

### 18S
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.tx.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132.wang.taxonomy)

get.current()
```

## Compute distance matrix and cluster OTUs
Calculate uncorrected pairwise distances between the aligned DNA sequences; distances > 0.1 are ignored to safe resources (note: avoid `cluster.split`, results are inaccurate). Then, cluster this sequences into OTUs with the [OptiClust](https://msphere.asm.org/content/2/2/e00073-17) method.
```bash
dist.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.1)

get.current()

cluster(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03, method=opti)

get.current()
```

## Get OTU table
Create a .shared and a .rabund file samplewise. Set the threshold (by convention a 3 % barcoding gap) and create the OTU-table.
```bash
make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)
```

## Count &alpha;-diversity
Count OTUs sample-wise with [AWK](https://en.wikipedia.org/wiki/AWK) (awk -v: var=val, -e: use program-text, OFS: output field separator, NF: input field number, NR: total number of input records so far).
```bash
system(awk -v OFS='\t' -e '{notus=0; for (i=4; i<=NF; i++) { if ($i > 0) notus++; }; if (NR > 1) print $2 OFS notus; }' wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.sample.summary)

get.current()
```

## Get representative sequences for OTUs
* = cluster centroids?

Returns a FASTA file where the headers additionally carry the OTU number and the total abundance.
```bash
get.oturep(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
```


## Classify clustered OTUs
Get consenus taxonomy for OTUs.
```bash
### 16S
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP16S.wang.taxonomy, label=0.03)

### 18S
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v132_EMP18S.wang.taxonomy, label=0.03)

get.current()
```

## Rarefaction curves
Calculate sample-wise rarefaction curves to control for sufficient sampling (here, sequencing) depth.
```bash
rarefaction.single(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, calc=sobs, freq=100)
```
Further process them with [this R script](mothur_RarefactionCurves.R).

![Rarefaction curves](/Graphs/mothur_RarefactionCurves.png)

## Track reads
Survey where the reads are 'lost' in the pipeline.
```bash
awk -f transposeList2Table.awk <(grep ’mothur > \|# of’ mothur.*.logfile | grep ’# of\|make.contigs\|screen.seqs\|align.seqs\|filter.seqs\|pre.cluster\|chimera.vsearch\|awk.\+single.accnos’ | sed -r ’s/mothur > /mothur:/’ | sed -r ’s/mothur/\nmothur/’ | sed -r ’s/\t/ /’ | sed -r ’s/#/Number/’) > countReadsMothur.csv
```
Further process them with [this R script](mothur_TrackReads.R).

![Reads loss in the mothur pipeline](/Graphs/mothur_TrackReads.png)
