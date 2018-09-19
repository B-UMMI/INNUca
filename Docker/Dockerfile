FROM ubuntu:16.04
MAINTAINER Miguel Machado <mpmachado@medicina.ulisboa.pt>
LABEL version="4.0-dev"

WORKDIR /NGStools/

RUN apt-get update

# -- General Dependencies ---
RUN apt-get install -y git wget make g++

# -- INNUca General Dependencies ---
RUN apt-get install -y python-dev default-jre

# --- kraken2 --
# -- dependencies --
RUN apt-get install -y libgomp1
# -- kraken2 --
RUN wget https://github.com/DerrickWood/kraken2/archive/v2.0.7-beta.tar.gz && \
    tar xf v2.0.7-beta.tar.gz && \
    cd kraken2-2.0.7-beta/ && \
    ./install_kraken2.sh /NGStools/kraken2-2.0.7-beta/ && \
    cd /NGStools && \
    rm -r v2.0.7-beta.tar.gz
ENV PATH="/NGStools/kraken2-2.0.7-beta:${PATH}"

# -- mlst Dependencies --
# Blast
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && \
    tar -xf ncbi-blast-2.6.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.6.0+-x64-linux.tar.gz
ENV PATH="/NGStools/ncbi-blast-2.6.0+/bin:${PATH}"
# Perl libs
RUN apt-get install -y libmoo-perl liblist-moreutils-perl libjson-perl libfile-slurp-perl

# --- mlst ----
RUN git clone https://github.com/tseemann/mlst.git
ENV PATH="/NGStools/mlst/bin:${PATH}"

# --- ReMatCh ---
RUN git clone https://github.com/B-UMMI/ReMatCh.git
ENV PATH="/NGStools/ReMatCh:${PATH}"

# --- INNUca ---
RUN git clone https://github.com/B-UMMI/INNUca.git
ENV PATH="/NGStools/INNUca:${PATH}"

# fixing permissions for MLST update
RUN chmod +x /NGStools/INNUca/Docker/update_mlst_db.sh && chmod o+wr /NGStools/mlst/scripts/ && chmod -R o+wr /NGStools/mlst/db/ && sed -i "s#OUTDIR=pubmlst#OUTDIR=/NGStools/mlst/scripts/pubmlst#1" /NGStools/mlst/scripts/mlst-download_pub_mlst

# Clean
RUN apt-get remove -y make g++ && \
    apt-get autoclean -y

WORKDIR /data/
