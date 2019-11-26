FROM perl:5.30-slim-stretch
MAINTAINER Miguel Machado <mpmachado@medicina.ulisboa.pt>
LABEL version="4.2.1-01"

WORKDIR /NGStools/

RUN apt-get update

# -- General Dependencies ---
RUN apt-get install -y git wget g++

# -- INNUca General Dependencies ---
RUN apt-get install -y python-dev python-pip python3 python3-pip procps libfontconfig1
# - Java -
RUN wget https://download.java.net/openjdk/jdk8u40/ri/openjdk-8u40-b25-linux-x64-10_feb_2015.tar.gz && \
    tar xf openjdk-8u40-b25-linux-x64-10_feb_2015.tar.gz && \
    rm openjdk-8u40-b25-linux-x64-10_feb_2015.tar.gz
ENV PATH="/NGStools/java-se-8u40-ri/bin:${PATH}"

RUN pip install plotly

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
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
    tar -xf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.9.0+-x64-linux.tar.gz
# any2fasta
RUN git clone https://github.com/tseemann/any2fasta.git
# Perl libs
RUN cpan Moo List::MoreUtils JSON

ENV PATH="/NGStools/ncbi-blast-2.9.0+/bin:/NGStools/any2fasta:${PATH}"

# --- mlst ----
RUN git clone https://github.com/tseemann/mlst.git
ENV PATH="/NGStools/mlst/bin:${PATH}"

# --- ReMatCh ---
# TODO: to be used after converting INNUca do Python v3
#RUN git clone https://github.com/B-UMMI/ReMatCh.git && \
#    cd ReMatCh && \
#    pip3 install setuptools && \
#    python3 setup.py install && \
#    cd /NGStools
RUN git clone https://github.com/B-UMMI/ReMatCh.git
ENV PATH="/NGStools/ReMatCh/ReMatCh/src/samtools-1.3.1/bin:/NGStools/ReMatCh/ReMatCh/src/bcftools-1.3.1/bin:/NGStools/ReMatCh/ReMatCh/src/bowtie2-2.2.9:/NGStools/ReMatCh/ReMatCh:${PATH}"

# --- INNUca ---
RUN git clone https://github.com/B-UMMI/INNUca.git && \
    pip install setuptools
ENV PATH="/NGStools/INNUca/src/fastqc_v0.11.5:/NGStools/INNUca/src/pilon_v1.23:/NGStools/INNUca/src/SPAdes-3.13.0-Linux/bin:/NGStools/INNUca/src/Trimmomatic-0.38:/NGStools/INNUca:${PATH}"

# fixing permissions for MLST update
RUN chmod +x /NGStools/INNUca/Docker/update_mlst_db.sh && chmod o+wr /NGStools/mlst/scripts/ && chmod -R o+wr /NGStools/mlst/db/ && sed -i "s#OUTDIR=pubmlst#OUTDIR=/NGStools/mlst/scripts/pubmlst#1" /NGStools/mlst/scripts/mlst-download_pub_mlst

# Clean
RUN pip3 uninstall setuptools && \
    apt-get remove -y g++ python-pip python3-pip && \
    apt-get autoclean -y

WORKDIR /data/
