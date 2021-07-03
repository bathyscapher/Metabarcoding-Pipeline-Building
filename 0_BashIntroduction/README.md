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

### FASTQ files
The data are hosted on a private NextCloud server (access will be provided only during the workshop). `<share id>` is the string after the s/ and before the /download in the link.
```
curl -u "shareid:password" -H 'X-Requested-With: XMLHttpRequest' 'https://nextcloud.somedomain.com/public.php/webdav/'
```

## Help
Each bash command has a help page:
```
man R
```






