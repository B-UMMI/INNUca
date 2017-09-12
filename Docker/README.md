INNUca.py - Docker
===============
INNUca - Reads Control and Assembly

*INNUENDO quality control of reads, de novo assembly and contigs quality assessment, and possible contamination search*

<https://github.com/B-UMMI/INNUca>


This is a dockerfile for using INNUca, with all dependencies already installed.

Within this container you can find:
- ubuntu:16.04
- git
- Python v2.7
- Java-JRE
- [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) v2.6.0
- [mlst](https://github.com/tseemann/mlst)
- [ReMatCh](https://github.com/B-UMMI/ReMatCh) v3.2
- [INNUca](https://github.com/B-UMMI/INNUca) v3.1



### Using play-with-docker
[![Try in PWD](https://cdn.rawgit.com/play-with-docker/stacks/cff22438/assets/images/button.png)](http://labs.play-with-docker.com/)

Within [play-with-docker](http://labs.play-with-docker.com/) webpage click on **create session**. Then, another page
will open with a big counter on the upper left corner. Click on **+ add new instance** and a terminal like instance should be generated on the right. On
this terminal you can load this docker image as follows:

`docker pull ummidock/innuca:3.1`

#### Build this docker on your local machine

For this, docker needs to be installed on your machine. Instructions for this can be found [here](https://docs.docker.com/engine/installation/).

##### Using DockerHub (automated build image)

`docker pull ummidock/innuca:3.1`

##### Using GitHub (build docker image)

1) `git clone https://github.com/B-UMMI/INNUca.git`  
2) `docker build -t innuca:3.1 ./INNUca/Docker/`

### Run (using automated build image)
    docker run --rm -u $(id -u):$(id -g) -it -v /local/folder/fastq_data:/data/ ummidock/innuca:3.1 INNUca.py --speciesExpected "Streptococcus agalactiae" --genomeSizeExpectedMb 2.1 --inputDirectory /data/ --outdir /data/innuca_output/ --threads 8 --maxNumberContigs 100



### Updating the mlst database in docker instance

After you've built the docker image, you can still update the mlst database. For this, the `update_mlst_db.sh` script is provided. Simply run in after initiating the instance with:

`/NGStools/INNUca/Docker/update_mlst_db.sh`

For more information on this please consult the [provided information](https://github.com/tseemann/mlst#updating-the-database) in the [mlst page](https://github.com/tseemann/mlst).

Contact
-------
Miguel Machado <mpmachado@medicina.ulisboa.pt>  
Catarina Mendes
<cimendes@medicina.ulisboa.pt>
