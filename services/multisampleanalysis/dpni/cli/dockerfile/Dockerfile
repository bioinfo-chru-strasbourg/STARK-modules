##########################################################################
# Dockerfile Version:   0.1
# Software:             DPNI
# Software Version:     0.1
# Description:          Dockerfile to build DPNI analysis container image
##########################################################################

ARG STARK_IMAGE="stark/stark:0.9.18.3"

########
# FROM #
########

FROM $STARK_IMAGE
LABEL maintainer="LAMOUCHE Jean-Baptiste <jean-baptiste.lamouche@chru-strasbourg.fr>" \
    Software="DPNI" \
    Version="0.1" \
    Description="DPNI"


ENV SRC=/src
RUN yum -y install  epel-release
RUN yum -y install  gcc gcc-c++ \
    bzip2 bzip2-devel \
    zlib2 zlib2-devel \
    xz xz-devel \
    unzip wget make docker


#############
# Snakemake #
#############
RUN cd /tmp  && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/lib/miniconda3/ && \
    echo 'export PATH="$PATH:/usr/local/lib/miniconda3/bin"' > /etc/profile.d/miniconda.sh && \
    chmod 755 /etc/profile.d/miniconda.sh && \
    export PATH="$PATH:/usr/local/lib/miniconda3/bin" && \
    /usr/local/lib/miniconda3/bin/conda install -y -c bioconda -c conda-forge snakemake bedtools samtools && \
    chmod -R 755 /usr/local/lib/miniconda3/ && \
    rm -f /tmp/Miniconda3-latest-Linux-x86_64.sh;

##########
#  java  #
##########
RUN yum -y install java java-1.8.0

##############
#  Bedtools  #
##############

ENV BEDTOOLS=/STARK/tools/bedtools/current/bin


########
# POOL #
########
COPY . /app

##########
### DPNI #
##########
#
RUN git clone https://github.com/JbaptisteLam/DPNI.git /app/lib/DPNI
ENV PATH="/usr/local/lib/miniconda3/bin:$PATH"
RUN conda env create -f /app/lib/DPNI/config/environment.yml
RUN echo "source activate dpni" >> ~/.bashrc

RUN chmod -R +x /app
ENV PATH="/app/bin:$PATH"

WORKDIR /STARK
ENTRYPOINT []

CMD ["/app/bin/pool","--help"]
