# Dockerfile
FROM centos:7

RUN yum -y install python27 git wget gcc gcc-c++ zlib-devel --nogpgcheck  

# PHP Extensions
RUN mkdir /gemini && \
    cd /gemini && \
    wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py && \
    python gemini_install.py /gemini /gemini/data --nodata

ENV PATH=$PATH:/gemini/bin

