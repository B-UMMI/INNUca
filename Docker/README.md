INNUca.py - Docker
===============
INNUca - Reads Control and Assembly

*INNUENDO quality control of reads, de novo assembly and contigs quality assessment, and possible contamination search*

<https://github.com/B-UMMI/INNUca>

## Dockerfile for INNUca

This is a dockerfile for using INNUca, with all dependencies already installed.

#### Using play-with-docker
[![Try in PWD](https://cdn.rawgit.com/play-with-docker/stacks/cff22438/assets/images/button.png)](http://labs.play-with-docker.com/)

Within [play-with-docker](http://labs.play-with-docker.com/) webpage click on **create session**. Then, another page
will open with a big counter on the upper left corner. Click on **+ add new instance** and a terminal like instance should be generated on the right. On
this terminal you can load this docker image as follows:

1) `docker pull cimendes/innuca`
2) `docker run -it cimendes/innuca bash`

#### Build this docker on your local machine

For this, docker needs to be installed on your machine. Instructions for this can be found [here](https://docs.docker.com/engine/installation/).

##### Using DockerHub

1) `docker pull cimendes/innuca`

To run, do
2) `docker run -it cimendes/innuca bash`

##### Using GitHub

1) `git clone https://github.com/cimendes/INNUca.git`
2) `docker build . INNUca`

To run, do
3) `docker run -it INNUca bash`


Contact
-------
Catarina Mendes
<cimendes@medicina.ulisboa.pt>



> Written with [StackEdit](https://stackedit.io/).
