FROM centos:7
MAINTAINER Victor Grentzinger <victor.grentzinger@gmail.com>

##################################################################################
# Dockerfile Version:   1.0
# Software:             listener
# Software Version:     1.0
# Description:          Dockerfile to build Microservices Listener container image
##################################################################################

ENV SRC=/src
RUN yum -y install  epel-release
RUN yum -y install  gcc gcc-c++ \
					bzip2 bzip2-devel \
					zlib2 zlib2-devel \
					xz xz-devel \
					unzip wget make \
					git libjpeg-turbo-devel docker


#############
#  python3  #
#############
RUN yum -y install https://repo.ius.io/ius-release-el7.rpm
RUN yum -y install python36u python36u-libs python36u-devel python36u-pip
RUN pip3.6 install snakemake

##########
#  java  #
##########
RUN yum -y install java java-1.8.0

#######
#  R  #
#######
RUN yum -y install R;
RUN Rscript -e 'install.packages("nnls",repos="http://cran.us.r-project.org")' ; \
	Rscript -e 'install.packages("Hmisc",repos="http://cran.us.r-project.org")' ; \
	Rscript -e 'install.packages("mgcv",repos="http://cran.us.r-project.org")' ; \
	Rscript -e 'install.packages("dplyr",repos="http://cran.us.r-project.org")'; \
	Rscript -e 'install.packages("plyr",repos="http://cran.us.r-project.org")';

#########
#  TCL  #
#########
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

##############
#  Bedtools  #
##############
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz ; \
	tar xf bedtools-2.27.1.tar.gz ; \
	rm -rf bedtools-2.27.1.tar.gz ; \
	cd bedtools2 ; \
	make prefix=$SRC/bedtools/2.27.1 install ; \
	cd ../ ; \
	rm -rf bedtools2;
ENV BEDTOOLS=$SRC/bedtools/2.27.1/bin

##########
#  GATK  #
##########
RUN  wget "https://github.com/broadinstitute/gatk/releases/download/4.1.4.0/gatk-4.1.4.0.zip" -O GenomeAnalysisTK-4.1.4.0.zip; \
	unzip GenomeAnalysisTK-4.1.4.0.zip ; \
	rm -rf GenomeAnalysisTK-4.1.4.0.zip ; \
	cd gatk-4.1.4.0 ; \
	mkdir -p $SRC/gatk/4.1.4.0/bin ; \
	cp gatk-package-4.1.4.0-local.jar $SRC/gatk/4.1.4.0/bin ; \
	cd ../ ; \
	rm -rf gatk-4.1.4.0;

#############
#  AnnotSV  #
#############
RUN wget https://lbgi.fr/AnnotSV/Sources/AnnotSV_2.1.tar.gz; \
		tar -xf AnnotSV_2.1.tar.gz; \
		rm -rf AnnotSV_2.1.tar.gz; \
		mkdir -p $SRC/AnnotSV/2.1; \
		rm -rf AnnotSV_2.1/Annotations; \
		mv AnnotSV_2.1/* $SRC/AnnotSV/2.1; \
		rm -rf AnnotSV_2.1;
ENV ANNOTSV=$SRC/AnnotSV/2.1

############
#  CANOES  #
############

COPY . /app
RUN echo -e '{\n\t"bedtools":"/src/bedtools/2.27.1/bin",\n\t"gatk":"/src/gatk/4.1.4.0/bin/gatk-package-4.1.4.0-local.jar",\n\t"java":"/usr/bin/java",\n\t"annotsv":"/src/AnnotSV/2.1",\n\t"R":"/usr/bin/R"\n}' > /app/config/config.json

RUN chmod -R +x /app
ENV PATH="/app/bin:$PATH"

WORKDIR /app

CMD ["bin/canoes","--help"]
