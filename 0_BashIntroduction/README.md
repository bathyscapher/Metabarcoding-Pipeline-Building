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

## 
```
```
