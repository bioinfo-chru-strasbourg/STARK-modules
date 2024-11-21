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
RUN dnf -y install wget rsync gcc-c++ which git python3.12 pandas

#bcftools
WORKDIR /tools
RUN wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
RUN tar -xvf bcftools-1.8.tar.bz2
RUN rm -rf bcftools-1.8.tar.bz2
WORKDIR /tools/bcftools-1.8
RUN ./configure --prefix=/tools/
RUN make
RUN make install

#tabix
WORKDIR /tools
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
RUN tar -xvf htslib-1.20.tar.bz2
RUN rm -rf htslib-1.20.tar.bz2
WORKDIR /tools/htslib-1.20
RUN ./configure --prefix=/tools/
RUN make
RUN make install

#bedtools
WORKDIR /tools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
RUN tar -xvf bedtools-2.29.1.tar.gz
RUN rm -rf bedtools-2.29.1.tar.gz
WORKDIR /tools/bedtools2
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN make

COPY bin/*.py /tools/vannot/
RUN chmod -R 777 /tools/

RUN echo 'alias ll="ls -lra"' >> ~/.bashrc
RUN echo 'alias vannot="python3 /tools/vannot/__main__.py"' >> ~/.bashrc
RUN echo 'PATH="$PATH:/tools/bin/"' >> ~/.bashrc
RUN echo 'PATH="$PATH:/tools/bedtools2/bin/"' >> ~/.bashrc

WORKDIR /tools/
RUN git clone https://github.com/SamuelNicaise/vannotplus
RUN /usr/bin/python3 -m ensurepip --upgrade
WORKDIR /tools/vannotplus
RUN /usr/bin/python3 -m pip install -e .
WORKDIR /