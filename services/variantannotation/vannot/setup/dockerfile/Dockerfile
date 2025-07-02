# source image
FROM centos:7
LABEL name="Mateusz RAUCH" email="<rauchmateusz@gmail.com>"

##########################################################################
# Dockerfile Version:   1.0
# Software:             vAnnot
# Software Version:     1.0
# Software Website:     http://www.lbgi.fr/VaRank/
# Description:          Dockerfile to build VaRank setup container image
##########################################################################

# packages installation
RUN yum -y update
RUN yum -y install epel-release
RUN yum -y install wget
RUN yum -y install make
# python, 
RUN yum -y install gcc
# python, 
RUN yum -y install bzip2 bzip2-devel
RUN yum -y install zlib-devel
# varank, 
RUN yum -y install git
# pandas
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
ENV PYTHONDONTWRITEBYTECODE=1

# ##########################################################################

# copy app binaries
COPY bin/setup.py /setup/setup.py
COPY config/* /setup/

# ##########################################################################

RUN chmod -R 777 /setup/

ENTRYPOINT ["python3.10 /setup/setup.py"]
WORKDIR /