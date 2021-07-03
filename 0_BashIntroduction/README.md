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
cd home
ll
cd rstudio
```

## Run a programs in the shell 
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
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
```


## Handle files
Bash supports any file handling.

1. First, navigate to the home directory
1. Create a directory
1. Enter this directory
1. Create a text file with the text editor [vim](https://github.com/vim/vim) and write some text into it (first, press 'i' to start editing mode, write the text, press `ESC` to return to normal mode, type `:x!` to safe and exit vim.
1. Use `more` or...
1. ... `head` to show the content of a file
1. Copy the file
1. To append text to a file use `echo` and `>>`
1. Show content of both files
1. Delete the first file
1. Navigate back to `home` ...
1. ... and delete the directory `test`

```
cd
mkdir test
cd test
vim testi.test
more testi.test
head testi.test
cp testi.test another.test
echo 'Some more text' >> testi.test
more *test
rm testi.test
cd
rm -rf test
```

Now, let's have a look at a bit more complex files.

1. Starting in `home`,
1. copy a *tar.gz file to home and ...
1. ... look into it
1. Extract and ...
1. list files
1. Look into the extracted files
1. Delete these files again

```
cd
cp ...tar.gz ~
more ...tar.gz
tar xzvf ...targ.z
ll
more ...*
rm ...
```


## Help
Each bash command has a help page:
```
man R
```




