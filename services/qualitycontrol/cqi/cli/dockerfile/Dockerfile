FROM continuumio/miniconda:latest

###############################
####  CQI Dockerfile
###############################

#necessary libs
RUN mkdir -p /usr/share/man/man1
RUN apt-get --allow-releaseinfo-change update && apt-get install -y \
	build-essential \
	default-jdk \
	git \
	vim \
	unzip \
	graphviz \
	gv \
	docker.io \
	bc \
	jq \
	rsync


#libs for htslib and  bcftools
RUN apt-get install -y zlib1g-dev 
RUN apt-get install -y liblzma-dev 
RUN apt-get install -y libbz2-dev
RUN apt-get install -y libcurl4-openssl-dev

#htslib
ENV HTSLIB_INSTALL_DIR=/opt/htslib
ENV HTSLIB_VERSION=1.14
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 && \
	tar --bzip2 -xf htslib-$HTSLIB_VERSION.tar.bz2
WORKDIR /tmp/htslib-$HTSLIB_VERSION
RUN make prefix=$HTSLIB_INSTALL_DIR && \
	make prefix=$HTSLIB_INSTALL_DIR install
WORKDIR /
RUN ln -s $HTSLIB_INSTALL_DIR/bin/bgzip /usr/bin/bgzip
RUN ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix
RUN rm -rf /tmp/htslib-$HTSLIB_VERSION

#bcftools
ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV BCFTOOLS_VERSION=1.14
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
	tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2
WORKDIR /tmp/bcftools-$BCFTOOLS_VERSION
RUN make prefix=$BCFTOOLS_INSTALL_DIR && \
	make prefix=$BCFTOOLS_INSTALL_DIR install
WORKDIR /
RUN ln -s $BCFTOOLS_INSTALL_DIR/bin/bcftools /usr/bin/bcftools && \
	rm -rf /tmp/bcftools-$BCFTOOLS_VERSION


#set up env 
WORKDIR /
COPY . /
COPY ./app/config/variant.yml /

#prepare env (base)
RUN conda env list

# create environment (variant)
RUN conda env create -f variant.yml

# activate environment
RUN echo "source activate variant" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

CMD ["/bin/bash"]

