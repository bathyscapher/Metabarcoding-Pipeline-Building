# Pipeline Building OTU
Generally following the [MiSeq SOP](https://mothur.org/wiki/miseq_sop/), the major steps are:
1. Read merging/contigs assembly
1. Alignment to reference database (SILVA)
1. Pre-clustering
1. Chimera detection
1. Singleton filtering
1. Distance matrix
1. OTU clustering and classification

Useful links are the [mothur wiki](https://mothur.org/wiki/), in particular the [command page](https://mothur.org/wiki/tags/#commands) and the [manual](https://mothur.org/wiki/mothur_manual/). Download the [most recent mothur version here](https://github.com/mothur/mothur/releases).


## mothur pipeline
### 16S & 18S
The pipeline is basically the same for 16S and 18S. Exceptions are steps dealing with the SILVA reference database that was customized to each primer pair and some consequential steps - such cases are highlighted. 

Open `mothur` directly in the shell (command line/terminal). Calling `mothur` opens the program and the commands below can be run from there.

Set working directory and specify number of available processors.
```bash
set.dir(tempdefault=.)
set.current(processors=6)
```
*Tip*: run `nproc` in the shell to see number of available processors.

### Merge reads to contigs
Create a metafile listing all the fastq in the directory with a given prefix (`wine.files`). Then, merge the forward and reverse reads to contigs with the [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).
```bash
make.file(inputdir=., type=gz, prefix=wine)
make.contigs(file=wine.files, align=needleman)

count.groups(group=wine.contigs.groups)
summary.seqs()
get.current()
```

### Quality trimming
Remove all contigs that are either too long or too short and exceed a certain amount of homopolymers. Values taken from output of previous `summary.seqs()`.

The expected amplicon length of the [EMP 16S primers](http://www.earthmicrobiome.org/protocols-and-standards/16s/) is about 300 to 350 bp, the [EMP 18S primers](http://www.earthmicrobiome.org/protocols-and-standards/18s/) about 260 +/- 50 bp.
```bash
### 16S
screen.seqs(fasta=wine.trim.contigs.fasta, group=wine.contigs.groups, summary=wine.trim.contigs.summary, maxambig=0, minlength=253, maxlength=254, maxhomop=8)

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

### Align the contigs to SILVA
Align the contigs to the (customized) [SILVA](https://www.arb-silva.de/) reference data base.

*Tip:* instructions to taylor the SILVA database to specific primers can be found [here](http://blog.mothur.org/2018/01/10/SILVA-v138.1-reference-files/) and [here](http://blog.mothur.org/2016/07/07/Customization-for-your-region/).

```bash
### 16S
align.seqs(fasta=wine.trim.contigs.good.unique.fasta, reference=../../silva/prok/silva.v138.1_16S-V4.align, flip=f)

### 18S
align.seqs(fasta=wine.trim.contigs.good.unique.fasta, reference=../../silva/euk/silva.v138.1_18S-V4.align, flip=f)

summary.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table)
get.current()
```

Screen once more to ensure that all sequences overlap the same region.
```bash
### 16S
screen.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table, summary=wine.trim.contigs.good.unique.summary, start=8, end=11570)

### 18S
screen.seqs(fasta=wine.trim.contigs.good.unique.align, count=wine.trim.contigs.good.count_table, summary=wine.trim.contigs.good.unique.summary, start=4, end=562)

summary.seqs(fasta=current, count=current)
```

Remove columns that contain only gaps (they stem from the alignment) and grab unique sequences only.
```bash
filter.seqs(fasta=wine.trim.contigs.good.unique.good.align, vertical=T, trump=.)
unique.seqs(fasta=wine.trim.contigs.good.unique.good.filter.fasta, count=wine.trim.contigs.good.good.count_table)

summary.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.fasta, count=wine.trim.contigs.good.unique.good.filter.count_table)
system(command='rm -f wine.filter')
```

### Precluster sequences
Unify similar sequences, i.e. cluster sequences that probably are noisy aka "denoising" (according to the MiSeq SOP, 1 difference per 100 bp is recommended).

**diff** | **2** | **3** | **5** | **6**
:--- | ---: | ---: | ---: | ---:
**nt for 2 seqs** | 4 | 6 | 10 | 12 
**% for 290 nt** | 1.4 | 2.1 | 3.4 | 4.1 

```bash
pre.cluster(fasta=wine.trim.contigs.good.unique.good.filter.unique.fasta, count=wine.trim.contigs.good.unique.good.filter.count_table, diffs=3, processors=6)

summary.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.count_table)
get.current()
```

### Detect chimeras
Detect chimeras with [vsearch](https://github.com/torognes/vsearch) and remove them.
```bash
chimera.vsearch(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
remove.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

count.groups(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
summary.seqs(fasta=current, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
get.current()
```

### Remove singletons
Generate a file containing all singletons with [AWK](https://en.wikipedia.org/wiki/AWK) and remove them from the current fasta and count_table.
```bash
system(awk -e '{if ($2 < 2) {print $1}}' wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table > wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)

remove.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)
remove.seqs(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, accnos=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.single.accnos)

count.groups(count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
summary.seqs(fasta=current, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
get.current()
```

### Classify sequences
Assign sequences to the reference database and taxonomy with the [Ribosomal Database Project (RDP) classifier](https://aem.asm.org/content/73/16/5261).
```bash
### 16S
classify.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, reference=../../silva/prok/silva.v138.1_16S-V4.align, taxonomy=../../silva/prok/silva.v138.1_16S-V4.tax, method=wang, cutoff=80)

### 18S
classify.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, reference=../../silva/prok/silva.v138.1_18S-V4.align, taxonomy=../../silva/prok/silva.v138.1_18S-V4.tax, method=wang, cutoff=80)

get.current()
```


### Compute distance matrix and cluster OTUs
Calculate uncorrected pairwise distances between the aligned DNA sequences; distances > 0.1 are ignored to safe resources (note: avoid `cluster.split`, results are inaccurate). Then, cluster this sequences into OTUs with the [OptiClust](https://msphere.asm.org/content/2/2/e00073-17) method.
```bash
dist.seqs(fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.1)
get.current()

cluster(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03, method=opti)
get.current()
```

### Get OTU table
Create a .shared and a .rabund file samplewise. Set the threshold (by convention a 3 % barcoding gap) and create the OTU-table.
```bash
make.shared(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
count.seqs(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)
```

### Classify clustered OTUs
Get consenus taxonomy for OTUs.
```bash
### 16S
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.1_16S_V4.wang.taxonomy, label=0.03)
         
### 18S
classify.otu(list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.v138.1_EMP18S.wang.taxonomy, label=0.03)

get.current()
```


### Get representative sequences for OTUs
Returns a FASTA file where the headers additionally carry the OTU number and the total abundance.
```bash
get.oturep(column=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, list=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, fasta=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=wine.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03)
```

### Rarefaction curves
Calculate sample-wise rarefaction curves to control for sufficient sampling (here, sequencing) depth.
```bash
rarefaction.single(shared=wine.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, calc=sobs, freq=100)
```
Further process them with [this R script](mothur_RarefactionCurves.R).

![Rarefaction curves](/Graphs/mothur_RarefactionCurves.png)


### Track reads
Survey where the reads are 'lost' in the pipeline.
```bash
awk -f /home/rstudio/Metabarcoding-Pipeline-Building/2_PipelineOTU/transposeList2Table.awk <(
grep 'mothur > \|# of' mothur.*.logfile | grep '# of\|make.contigs\|screen.seqs\|align.seqs\|filter.seqs\|pre.cluster\|chimera.vsearch\|awk.\+single.accnos' | sed -r 's/mothur > /mothur: /' | sed -r 's/mothur/\nmothur/' | sed -r 's/\t/ /' | sed -r 's/#/Number/') > countReads.csv
```
Further process them with [this R script](mothur_TrackReads.R).

![Read loss in the mothur pipeline](/Graphs/mothur_TrackReads.png)




