[![dockeri.co](https://dockeri.co/image/ummidock/innuca)](https://hub.docker.com/r/ummidock/innuca)

# INNUca.py - Docker

INNUca - Reads Control and Assembly

*INNUENDO quality control of reads, de novo assembly and contigs quality assessment, and possible contamination search*

<https://github.com/B-UMMI/INNUca>

<https://hub.docker.com/r/ummidock/innuca>


This is a dockerfile for using INNUca, with all dependencies already installed.

Within this container you can find:
- Debian Stretch (9)
- Perl v5.30
- git v2.11.0
- Python v2.7
- Java-JDK v1.8.0_40 headless
- [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) v2.9.0
- [mlst](https://github.com/tseemann/mlst) v2.18.0
- [ReMatCh](https://github.com/B-UMMI/ReMatCh) v4.1.0
- [Kraken](https://ccb.jhu.edu/software/kraken/) v2.0.7
- [INNUca](https://github.com/B-UMMI/INNUca) v4.2.2



### Using play-with-docker
[![Try in PWD](https://cdn.rawgit.com/play-with-docker/stacks/cff22438/assets/images/button.png)](http://labs.play-with-docker.com/)

Within [play-with-docker](http://labs.play-with-docker.com/) webpage click on **create session**. Then, another page
will open with a big counter on the upper left corner. Click on **+ add new instance** and a terminal like instance should be generated on the right. On
this terminal you can load this docker image as follows:

`docker pull ummidock/innuca:4.2.2-01`

#### Build this docker on your local machine

For this, docker needs to be installed on your machine. Instructions for this can be found [here](https://docs.docker.com/engine/installation/).

##### Using DockerHub (automated build image)

`docker pull ummidock/innuca:4.2.2-01`

##### Using GitHub (build docker image)

1) `git clone https://github.com/B-UMMI/INNUca.git`  
2) `docker build -t ummidock/innuca:4.2.2-01 ./INNUca/Docker/`

### Run (using automated build image)
    docker run --rm -u $(id -u):$(id -g) -it -v /local/folder/fastq_data:/data/ ummidock/innuca:4.2.2-01 INNUca.py --speciesExpected "Streptococcus agalactiae" --genomeSizeExpectedMb 2.1 --inputDirectory /data/ --outdir /data/innuca_output/ --threads 8 --maxNumberContigs 100

### udocker

> "A basic user tool to execute simple docker containers in user space without requiring root privileges.". From [here](https://github.com/indigo-dc/udocker).

```bash
# Get Docker image
udocker pull ummidock/innuca:4.2.2-01

# Create container (only needed to be done once)
udocker create --name=innuca_4-2-2_01 ummidock/innuca:4.2.2-01

# Run INNUca
udocker run --user $(id -u):$(id -g) -v /local/folder/fastq_data:/data/ innuca_4-2-2_01 INNUca.py --speciesExpected "Streptococcus agalactiae" --genomeSizeExpectedMb 2.1 --inputDirectory /data/ --outdir /data/innuca_output/ --threads 8 --maxNumberContigs 100
```
More examples on how to use **udocker** can be found in **udocker** [GitHub page](https://github.com/indigo-dc/udocker)  
  
*__NOTE__*: if some `Error: in download: HTTP/1.1 400 Bad Request` occur while pulling the Docker image, first pull it using `docker` and then retry pulling it with `udocker` (same with Shifter).

### Updating the mlst database in docker instance

After you've built the docker image, you can still update the mlst database. For this, the `update_mlst_db.sh` script is provided. Simply run in after initiating the instance with:

`/NGStools/INNUca/Docker/update_mlst_db.sh`

For more information on this please consult the [provided information](https://github.com/tseemann/mlst#updating-the-database) in the [mlst page](https://github.com/tseemann/mlst).

Contact
-------
Miguel Machado <mpmachado@medicina.ulisboa.pt>  
Catarina Mendes
<cimendes@medicina.ulisboa.pt>
