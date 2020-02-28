This repository creates a docker image that can be used to make dada2 analyses on fastq files.

# Dependencies

* Docker

# Using docker image

If you are just interested in running an analysis, you can download the docker image from docker hub:

`TODO: put it on docker hub`

The docker images expects to find the fastq input files in /media/input and will write the output to /media/output (inside your docker image). So mount the folder with the input files to /media/input and mount an empty folder to /media/output, like this:

```
docker run --mount type=bind,source=`pwd`/input-folder/,destination=/media/input --mount type=bind,source=`pwd`/output-folder/,destination=/media/output --mount type=bind,source=`pwd`/test-silva/,destination=/media/silva dada2-analysis --truncate_length_fwd 245 --truncate_length_rev 165 --truncate_quality 12 --rank Family --path_to_silva /media/silva/silva_nr_v128_train_set.fa.gz
```

The following arguments are available for tweaking the analysis:

* truncate\_length\_fwd
* truncate\_length\_rev
* truncate\_quality
* rank

Note: The mount commands need to be before the name of the docker image name, dada2-analysis, and all the other arguments needs to follow after the docker image name.

# Build from source

Build the docker image by running this command:

`docker build -t dada2-analysis .`

# Publish image

```
docker tag [SOURCE_IMAGE] [HOSTNAME]/[PROJECT-ID]/[IMAGE][:TAG]
gcloud docker -- push [HOSTNAME]/[PROJECT-ID]/[IMAGE][:TAG]
```

For example:

```
docker tag dada2-analysis:test eu.gcr.io/adroit-medium-188915/dada2-analysis:test
gcloud docker -- push eu.gcr.io/adroit-medium-188915/dada2-analysis:test
```
