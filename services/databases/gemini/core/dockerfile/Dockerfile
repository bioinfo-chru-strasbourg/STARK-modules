# Dockerfile
FROM stark/stark:0.9.18.2

RUN yum -y install python27 git wget gcc gcc-c++ zlib-devel --nogpgcheck  

# PHP Extensions
RUN mkdir -p /STARK/tools/gemini/latest && \
    cd /STARK/tools/gemini/latest && \
    wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py && \
    python gemini_install.py /STARK/tools/gemini/latest /STARK/tools/gemini/latest/data --nodata

ENV PATH=$PATH:/gemini/bin

