# source image
FROM almalinux:8
LABEL name="Mateusz RAUCH" email="<rauchmateusz@gmail.com>"

##########################################################################
# Dockerfile Version:   1.0
# Software:             VariantAnnotation
# Software Version:     1.0
# Software Website:     https://github.com/bioinfo-chru-strasbourg/howard
# Description:          Dockerfile to build VariantAnnotation image
##########################################################################

# DNF envs
RUN dnf -y install dnf-plugins-core
RUN dnf config-manager --add-repo=https://download.docker.com/linux/centos/docker-ce.repo
RUN dnf -y install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
RUN dnf -y install wget

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

COPY bin/*.py /tools/variantannotation/
RUN chmod -R 777 /tools/

RUN echo "source activate common" >> ~/.bashrc
RUN echo "alias ll='ls -lra'" >> ~/.bashrc
# RUN echo "alias variantannotation='python3.12 /tools/variantannotation/current/scripts/__main__.py'" >> ~/.bashrc

WORKDIR /