FROM ubuntu:16.04

WORKDIR /NGStools/

# -- General Dependencies ---
RUN apt-get update
RUN apt-get install -y git
RUN apt-get install -y python-dev
RUN apt-get install -y wget
RUN apt install -y openjdk-8-jre

# -- mlst Dependencies --
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
RUN tar -xf ncbi-blast-2.6.0+-x64-linux.tar.gz
RUN rm ncbi-blast-2.6.0+-x64-linux.tar.gz
ENV PATH="/NGStools/ncbi-blast-2.6.0+/bin:${PATH}"
RUN apt-get install -y libmoo-perl liblist-moreutils-perl

# --- mlst ----
RUN git clone https://github.com/tseemann/mlst.git
ENV PATH="/NGStools/mlst/bin:${PATH}"

# --- INNUca ---
#RUN git clone https://github.com/B-UMMI/INNUca.git
RUN git clone https://github.com/cimendes/INNUca.git

ENV PATH="/NGStools/INNUca:${PATH}"

#fixing permissions
RUN chmod +x /NGStools/INNUca/Docker/update_mlst_db.sh && chmod o+wr /NGStools/mlst/scripts/ && chmod -R o+wr /NGStools/mlst/db/ && sed -i "s#OUTDIR=pubmlst#OUTDIR=/NGStools/mlst/scripts/pubmlst#1" /NGStools/mlst/scripts/mlst-download_pub_mlst

#RUN INNUca.py --help

