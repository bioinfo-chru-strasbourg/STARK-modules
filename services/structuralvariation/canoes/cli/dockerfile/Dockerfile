##########################################################################
# Dockerfile Version:   1.1.1
# Software:             CANOES
# Software Version:     1.1.1
# Software Website:     http://www.columbia.edu/~ys2411/canoes/CANOES.R
# Description:          Dockerfile to build CANOES container image
##########################################################################


########
# FROM #
########

FROM centos:7
LABEL maintainer="LAMOUCHE Jean-Baptiste <jbaptiste.lamouche@gmail.com>" \
	Software="CANOES" \
	Version="1.1.1" \
	Description="CANOES"


ENV SRC=/tools
RUN yum -y install  epel-release  && \
	yum -y install  gcc gcc-c++ \
	bzip2 bzip2-devel \
	zlib2 zlib2-devel \
	xz xz-devel \
	unzip wget make bgzip tabix\
	ca-certificates \
	git libjpeg-turbo-devel docker bcftools \
	htslib-tools


#############
#  python3  #
#############
#RUN yum -y install https://repo.ius.io/ius-release-el7.rpm  && \
#	yum -y install python36u python36u-libs python36u-devel python36u-pip

#############
# Snakemake #
#############
#RUN cd /tmp  && \
#	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  && \
#	bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/lib/miniconda3/ && \
#	echo 'export PATH="$PATH:/usr/local/lib/miniconda3/bin"' > /etc/profile.d/miniconda.sh && \
#	chmod 755 /etc/profile.d/miniconda.sh && \
#	export PATH="$PATH:/usr/local/lib/miniconda3/bin" && \
#	/usr/local/lib/miniconda3/bin/conda install -y -c bioconda -c conda-forge snakemake && \
#	chmod -R 755 /usr/local/lib/miniconda3/ && \
#	rm -f /tmp/Miniconda3-latest-Linux-x86_64.sh;

##########
#java & R#
##########
RUN yum -y install java java-1.8.0 R
RUN Rscript -e 'install.packages("nnls",repos="http://cran.us.r-project.org")' && \
	Rscript -e 'install.packages("Hmisc",repos="http://cran.us.r-project.org")' && \
	Rscript -e 'install.packages("mgcv",repos="http://cran.us.r-project.org")' && \
	Rscript -e 'install.packages("dplyr",repos="http://cran.us.r-project.org")'&& \
	Rscript -e 'install.packages("plyr",repos="http://cran.us.r-project.org")'&& \
	Rscript -e 'install.packages("ggplot2",repos="http://cran.us.r-project.org")'&& \
	Rscript -e 'install.packages("reshape2",repos="http://cran.us.r-project.org")'&& \
	Rscript -e 'install.packages("gtools",repos="http://cran.us.r-project.org")'&& \
	Rscript -e 'install.packages("stringr",repos="http://cran.us.r-project.org")';

#########
#  TCL  #
#########
RUN wget --no-check-certificate "https://sourceforge.net/projects/tcl/files/Tcl/8.6.9/tcl8.6.9-src.tar.gz/download?use_mirror=netcologne" && \
	tar -xf download?use_mirror=netcologne && \
	rm -f download?use_mirror=netcologne && \
	cd tcl8.6.9/unix && \
	./configure --prefix=$SRC/tcl && \
	make && \
	make prefix=$SRC/tcl install && \
	ln -s tclsh8.6 $SRC/tcl/bin/tclsh && \
	cd ../../ && \
	rm -rf tcl8.6.9;
ENV PATH="$SRC/tcl/bin:$PATH"

##############
#  Bedtools  #
##############
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz  && \
	tar xf bedtools-2.27.1.tar.gz  && \
	rm -rf bedtools-2.27.1.tar.gz  && \
	cd bedtools2  && \
	make prefix=$SRC/bedtools/2.27.1 install  && \
	cd ../  && \
	rm -rf bedtools2;
ENV BEDTOOLS=$SRC/bedtools/2.27.1/bin

##########
#  GATK  #
##########
RUN  wget "https://github.com/broadinstitute/gatk/releases/download/4.1.4.0/gatk-4.1.4.0.zip" -O GenomeAnalysisTK-4.1.4.0.zip && \
	unzip GenomeAnalysisTK-4.1.4.0.zip  && \
	rm -rf GenomeAnalysisTK-4.1.4.0.zip  && \
	cd gatk-4.1.4.0  && \
	mkdir -p $SRC/gatk/4.1.4.0/bin  && \
	cp gatk-package-4.1.4.0-local.jar $SRC/gatk/4.1.4.0/bin  && \
	cd ../  && \
	rm -rf gatk-4.1.4.0;

#############
#  AnnotSV  #
#############
RUN wget https://lbgi.fr/AnnotSV/Sources/AnnotSV_2.1.tar.gz && \
	tar -xf AnnotSV_2.1.tar.gz && \
	rm -rf AnnotSV_2.1.tar.gz && \
	mkdir -p $SRC/AnnotSV/2.1 && \
	rm -rf AnnotSV_2.1/Annotations && \
	mv AnnotSV_2.1/* $SRC/AnnotSV/2.1 && \
	rm -rf AnnotSV_2.1 \
	rm -rf $SRC/AnnotSV/2.1/Annotations/Genes-based/OMIM/*
RUN mkdir -p $SRC/AnnotSV/2.1

#ENV ANNOTSV=$SRC/AnnotSV/2.1

###################
#  AnnotSV 3.1 #
###################

RUN mkdir -p $SRC/AnnotSV/3.1
WORKDIR $SRC/AnnotSV/3.1
RUN git clone -b v3.1 https://github.com/lgmgeo/AnnotSV.git $SRC/AnnotSV/3.1
RUN make PREFIX=. install
ENV ANNOTSV=$SRC/AnnotSV/3.1

###########################
## Miniconda & Snakemake  #
###########################

RUN cd /tmp  && \
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	sh Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/lib/miniconda3
ENV PATH="/usr/local/lib/miniconda3/bin:$PATH"
RUN source ~/.bashrc && conda install -y -c bioconda -c conda-forge snakemake

############
#  CANOES  #
############

COPY . /app
#RUN echo -e '{\n\t"bedtools":"/src/bedtools/2.27.1/bin",\n\t"gatk":"/src/gatk/4.1.4.0/bin/gatk-package-4.1.4.0-local.jar",\n\t"java":"/usr/bin/java",\n\t"annotsv":"/src/AnnotSV/2.1",\n\t"R":"/usr/bin/R"\n}' > /app/config/config.json
RUN conda env create -f /app/config/ngs.yaml
RUN echo "source activate ngs" >> ~/.bashrc
ENV PATH="/opt/conda/envs/env/bin:$PATH"
RUN chmod -R +x /app
ENV PATH="/app/bin:$PATH"

WORKDIR /app

############
#  CLEAN   #
############
RUN yum clean all 

CMD ["cat", "/app/config/default.yaml"]
