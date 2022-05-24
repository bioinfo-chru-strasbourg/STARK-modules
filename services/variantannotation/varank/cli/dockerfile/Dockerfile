# source image
FROM centos:7
LABEL name="Mateusz RAUCH" email="<rauchmateusz@gmail.com>"

##########################################################################
# Dockerfile Version:   1.0
# Software:             VaRank
# Software Version:     2.0
# Software Website:     http://www.lbgi.fr/VaRank/
# Description:          Dockerfile to build VaRank container image
##########################################################################

# packages installation
RUN yum -y update
RUN yum -y install epel-release
RUN yum -y install wget
RUN yum -y install make
# python
RUN yum -y install gcc
# python
RUN yum -y install bzip2 bzip2-devel
RUN yum -y install zlib-devel
RUN yum -y install libffi-devel
# varank
RUN yum -y install git
# panda
RUN yum clean all

# ##########################################################################

# OPENSSL
WORKDIR /usr/src
RUN wget --no-check-certificate "https://www.openssl.org/source/openssl-1.1.1o.tar.gz"
RUN tar xzf openssl-1.1.1o.tar.gz; \
	rm -f openssl-1.1.1o.tar.gz
WORKDIR /usr/src/openssl-1.1.1o
RUN ./config \
	--prefix=/usr/local/custom-openssl \
	--libdir=lib \
	--openssldir=/etc/pki/tls
RUN make -j1 depend
RUN make -j8
RUN make install_sw

# PYTHON
WORKDIR /usr/src
RUN wget --no-check-certificate "https://www.python.org/ftp/python/3.10.4/Python-3.10.4.tgz"
RUN	tar xzf Python-3.10.4.tgz; \
	rm -f Python-3.10.4.tgz
WORKDIR /usr/src/Python-3.10.4
RUN	./configure -C \
	--with-openssl=/usr/local/custom-openssl \
	--with-openssl-rpath=auto \
	--prefix=/usr/local/python-3.10.4
RUN make -j8
RUN	make altinstall
ENV PATH="/usr/local/python-3.10.4/bin:$PATH"

RUN pip3.10 install pandas
RUN pip3.10 install requests
ENV PYTHONDONTWRITEBYTECODE=1



# TCL
WORKDIR /usr/src
RUN wget --no-check-certificate "https://sourceforge.net/projects/tcl/files/Tcl/8.6.12/tcl8.6.12-src.tar.gz/download?use_mirror=netcologne"
RUN tar -xf download?use_mirror=netcologne; \
	rm -f download?use_mirror=netcologne
WORKDIR /usr/src/tcl8.6.12/unix
RUN ./configure --prefix=/usr/src/tcl; \
	make; \
	make prefix=/usr/src/tcl install;
RUN	ln -s tclsh8.6.12 /usr/src/tcl/bin/tclsh; \
	cd ../../; \
	rm -rf tcl8.6.12
ENV PATH="/usr/src/tcl/bin:$PATH"

# Java
RUN yum -y install java-11-openjdk

# SQLite3

RUN yum -y install sqlite-devel

# VaRank
RUN mkdir /tools; \
	mkdir /tools/varank
WORKDIR /tools/varank
RUN git clone https://github.com/lgmgeo/VaRank.git
RUN	ln -s VaRank current

ENV VARANK=/tools/varank/current/
ENV varank="python3 /tools/varank/current/scripts/__main__.py"

# Alamut-Batch
RUN mkdir /tools/alamut_batch
WORKDIR /tools/alamut_batch
RUN wget --no-check-certificate "https://downloads.interactive-biosoftware.com/download.php?product=Alamut-batch-standalone&platform=&OS=Linux&version=1.11"
RUN tar -xvf download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11; \
	rm -rf download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11; \
	ln -s alamut-batch-standalone-1.11 current
WORKDIR /tools/alamut_batch/current
RUN mv ancillary ancillary.original

COPY src/ancillary.tar.gz /tools/alamut_batch/current/
RUN tar -xvf ancillary.tar.gz
RUN chmod -R 777 /tools/alamut_batch/current/ancillary/
RUN rm -rf alamut-batch.ini
ENV ALAMUT=/tools/alamut_batch/current/

# ##########################################################################

# copy app binaries
COPY bin/*.py /tools/varank/current/scripts/

# ##########################################################################

RUN chmod -R 777 /tools/

ENTRYPOINT ["python3.10 /tools/varank/current/scripts/__main__.py"]
WORKDIR /