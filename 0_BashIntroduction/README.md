# Introduction to Bash

## Open Bash in Docker
```
docker container ls -a
docker exec -ti <container ID> bash
```

## Explore the directories
1. Show the present working directory (it is the "root" directory, symbolized by `/`):
1. List content in the current directory (works also with `ls`)
1. Change to the `home` directory (usually the space for the user in Linux)
1. List content once more: there are three directories:
    * `.` means the current directory
    * `..` means the parent directory
    * `rstudio` is the name of a directory
1. Change to `rstudio` ...
1. ... and list content again

```
pwd
ll
cd /home
ll
cd rstudio
ll
```

## Run a program in the shell 
### Open R
R runs directly in the shell:
```
R
```

### Show version
```
R --version
```

### Show path
```
which R
```


## Git
[Git](https://git-scm.com/) is a version control system that is employed by [GitHub](https://github.com/). The course material is all hosted on in a [GitHub repo(sitory)](https://github.com/bathyscapher/). Let's fetch this repo.
```
git clone https://github.com/bathyscapher/Metabarcoding-Pipeline-Building
```


## Downloads
### Reference database
Download the DADA2 reference database data directly from the shell into a target directory (`mothur/refs`). Then, check out this directory.
```
wget -P /mothur/refs https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz
cd /mothur/refs
ll -rt
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
```


## Handle files
Now, let's handle some file in Bash.

1. First, navigate to the home directory
1. Create a directory
1. Enter this directory
1. Echo some text into a not-yet-existing file.
1. Use `more` or...
1. ... `head` to show the content of a file
1. Copy the file
1. To append text to a file use `echo` and `>>`
1. Show content of both files
1. Delete the first file
1. Navigate back to `home` ...
1. ... and delete the directory `test`

```
cd /home/
mkdir test
cd test
echo 'This is any text.' > testi.test
more testi.test
head testi.test
cp testi.test another.test
echo 'Some more text' >> testi.test
more *test
rm testi.test
cd ..
rm -r test
```

Now, let's have a look at a bit more complex files.

1. Starting in `home`,
1. copy a *tar.gz file to home (the `.` stands for the present working directory) and ...
1. ... look into it ...
1. ... and exit with `q`.
1. Extract and ...
1. ...list files with the `-h` flag (for human readable).
1. `touch` one of the files and ...
1. ... listing the files shows an updated time stamp.
1. Look into the extracted file ...
1. ... and quit with `q`.
1. To only see the headers of the FASTA, use `grep` for extracting all lines with a '>'
1. Delete the file again.

```
cd /home/
cp /mothur/refs/silva_nr_v132_train_set.fa.gz .
more silva_nr_v132_train_set.fa.gz
q
gunzip silva_nr_v132_train_set.fa.gz
ll -h
touch silva_nr_v132_train_set
ll
more silva_nr_v132_train_set
q
grep '>' silva_nr_v132_train_set.fa
grep Negativicutes silva_nr_v132_train_set.fa
rm silva_nr_v132_train_set
```


## Help
Each bash command has a help page:
```
man R
```




