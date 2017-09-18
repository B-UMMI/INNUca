FROM ubuntu:16.04

WORKDIR /NGStools/

RUN apt-get update

# -- General Dependencies ---
RUN apt-get install -y git wget

# -- INNUca General Dependencies ---
RUN apt-get install -y python-dev default-jre

# -- mlst Dependencies --
# Blast
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && tar -xf ncbi-blast-2.6.0+-x64-linux.tar.gz && rm ncbi-blast-2.6.0+-x64-linux.tar.gz
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
RUN apt-get remove -y wget && apt-get autoclean -y

WORKDIR /data/
