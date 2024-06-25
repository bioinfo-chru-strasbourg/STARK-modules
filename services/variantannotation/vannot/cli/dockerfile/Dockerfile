# source image
FROM almalinux:8
LABEL name="Mateusz RAUCH" email="<rauchmateusz@gmail.com>"

##########################################################################
# Dockerfile Version:   1.0
# Software:             vannot
# Software Version:     1.0
# Software Website:     https://github.com/bioinfo-chru-strasbourg/howard
# Description:          Dockerfile to build vannot image
##########################################################################

# DNF envs
RUN dnf -y install dnf-plugins-core
RUN dnf config-manager --add-repo=https://download.docker.com/linux/centos/docker-ce.repo
RUN dnf -y install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
RUN dnf -y install gcc zlib-devel bzip2 bzip2-devel xz-devel make
RUN dnf -y install wget rsync

#mamba
WORKDIR /tmp
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
RUN chmod 777 Mambaforge-Linux-x86_64.sh
RUN ./Mambaforge-Linux-x86_64.sh -b -p /usr/local/lib/mamba/
RUN rm -rf Mambaforge-Linux-x86_64.sh
ENV MAMBA="/usr/local/lib/mamba/bin/mamba"

#work env
RUN $MAMBA create -n common python=3.12
RUN $MAMBA init

#bcftools
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
RUN tar -xvf bcftools-1.8.tar.bz2
WORKDIR /tmp/bcftools-1.8
RUN ./configure --prefix=/tools/
RUN make
RUN make install
ENV bcftools=/tools/bin/bcftools

#tabix
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
RUN tar -xvf htslib-1.20.tar.bz2
WORKDIR /tmp/htslib-1.20
RUN ./configure --prefix=/tools/
RUN make
RUN make install
ENV tabix=/tools/bin/tabix

COPY bin/*.py /tools/vannot/
RUN chmod -R 777 /tools/

RUN echo "source activate common" >> ~/.bashrc
RUN echo "alias ll='ls -lra'" >> ~/.bashrc
# RUN echo "alias vannot='python3.12 /tools/vannot/current/scripts/__main__.py'" >> ~/.bashrc

WORKDIR /