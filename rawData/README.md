Relevant Workflow for Getting Data
==================================

## Usage ##
Exverything is wrapped into a nice little make file

```
make getData
```

## Outline ##
First we need to find everything relevant in `/data/archive/` and copy it here

```
find /data/archive/ectemp/error\ correction/sequencing/miseq2/02* -type f -name "*.fastq.gz" | xargs -i cp {} .

```

Next, we will clean all the file names with a `bash` script.

```
./cleanNames.sh
```

This will get the files to a reasonable place, but since they are so mangled to begin with, we are going to have do do some editing by hand... :C Fortunately for you, I took care of this already!

```
./handClean.sh
```

Everything should be good to go!
