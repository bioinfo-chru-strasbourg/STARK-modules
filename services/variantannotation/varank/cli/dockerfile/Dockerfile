# source image
FROM centos:7
MAINTAINER Mateusz RAUCH <rauchmateusz@gmail.com>

##########################################################################
# Dockerfile Version:   1.0
# Software:             VaRank
# Software Version:     1.4.3
# Software Website:     http://www.lbgi.fr/VaRank/
# Description:          Dockerfile to build VaRank container image
##########################################################################

ENV SRC=/src
# packages installation
RUN yum -y update && yum -y install  epel-release && yum -y group install 'Development Tools' \
    && yum -y install  gcc gcc-c++ \
					bzip2 bzip2-devel \
					zlib2 zlib2-devel \
					xz xz-devel \
					unzip wget make \
					git libjpeg-turbo-devel \
					openssl-devel libffi-devel \
					python3.8.5-pip \
	&& yum clean all

##########################################################################

# tools installation
# PYTHON 3.8.5
RUN cd /usr/src; \
	wget "https://www.python.org/ftp/python/3.8.5/Python-3.8.5.tgz"; \
	tar xzf Python-3.8.5.tgz; \
	rm -f Python-3.8.5.tgz; \
	cd Python-3.8.5; \
	./configure --enable-optimizations; \
	make install;

RUN pip3 install pandas
				
# TCL
RUN wget "https://sourceforge.net/projects/tcl/files/Tcl/8.6.9/tcl8.6.9-src.tar.gz/download?use_mirror=netcologne"; \
	tar -xf download?use_mirror=netcologne; \
	rm -f download?use_mirror=netcologne; \
	cd tcl8.6.9/unix; \
	./configure --prefix=$SRC/tcl; \
	make; \
	make prefix=$SRC/tcl install; \
	ln -s tclsh8.6 $SRC/tcl/bin/tclsh; \
	cd ../../; \
	rm -rf tcl8.6.9;
ENV PATH="$SRC/tcl/bin:$PATH"

# Java
RUN yum -y install java-11-openjdk

# VaRank
RUN mkdir /tools; \
	mkdir /tools/varank; \
	cd /tools/varank; \
	wget "http://www.lbgi.fr/VaRank/Sources/VaRank_1.4.3.tar.gz"; \
	tar -xvf VaRank_1.4.3.tar.gz; \
	rm -rf VaRank_1.4.3.tar.gz; \
	ln -s VaRank_1.4.3 current; \
	cd current/ExtAnn; \
	rm -f OMIMannotations.tsv.gz; \
	ln -s /STARK/databases/OMIMannotations/current/OMIMannotations.tsv.gz ./OMIMannotations.tsv.gz

ARG VARANK=/tools/varank/current/
ENV VARANK $VARANK

# Alamut-Batch
RUN mkdir /tools/alamut_batch; \
	cd /tools/alamut_batch; \
	wget "https://downloads.interactive-biosoftware.com/download.php?product=Alamut-batch-standalone&platform=&OS=Linux&version=1.11"; \
	tar -xvf download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11; \
	rm -rf download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11; \
	ln -s alamut-batch-standalone-1.11 current; \
	cd current; \
	mv ancillary ancillary.original;

COPY src/ancillary /tools/alamut_batch/current/
ARG ALAMUT=/tools/alamut_batch/current/
ENV ALAMUT $ALAMUT

##########################################################################

# copy app binaries
COPY bin/newVaRank.py /tools/varank/current/scripts/newVaRank.py

##########################################################################

# copy config files

COPY config/* /tools/varank/current/config/

RUN chmod -R 777 /tools/

ENTRYPOINT ["python3", "/tools/varank/current/scripts/newVaRank.py"]
WORKDIR /