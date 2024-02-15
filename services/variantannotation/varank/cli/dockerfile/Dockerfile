# source image
FROM centos:7
LABEL name="Mateusz RAUCH" email="<rauchmateusz@gmail.com>"

##########################################################################
# Dockerfile Version:   2.1
# Software:             VaRank
# Software Version:     1.5.2
# Software Website:     http://www.lbgi.fr/VaRank/
# Description:          Dockerfile to build VaRank container image
##########################################################################

# packages installation
RUN yum -y update
RUN yum -y install epel-release
RUN yum -y install wget
RUN yum -y install htslib-tools
RUN yum -y install make
RUN yum -y install gcc
RUN yum -y install zlib-devel
RUN yum -y install bzip2-devel
RUN yum -y install xz-devel
# varank
RUN yum -y install git
RUN yum -y install glibc-2.17-326.el7_9.i686
RUN yum clean all

#conda
WORKDIR tmp
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/lib/miniconda3/
RUN echo 'export PATH="$PATH:/usr/local/lib/miniconda3/bin"' > /etc/profile.d/miniconda.sh
RUN chmod 755 /etc/profile.d/miniconda.sh
ENV PATH="$PATH:/usr/local/lib/miniconda3/bin"
RUN chmod -R 755 /usr/local/lib/miniconda3/
RUN rm -f /tmp/Miniconda3-latest-Linux-x86_64.sh

#work env
RUN conda install -c intel tcl
RUN conda create -n common python=3.12

# Java
RUN yum -y install java-11-openjdk

# VaRank
RUN mkdir /tools; \
	mkdir /tools/varank
WORKDIR /tools/varank
RUN git clone https://github.com/lgmgeo/VaRank.git
RUN ln -s VaRank current
RUN cd current; \
	git checkout c137f3a28c873d6e0cf66d71c84bd49765114031
RUN mkdir /tools/varank/current/tcllib ; \
	mkdir /tools/varank/current/Annotations_Exomiser
COPY src/tcllib.tar.gz /tools/varank/current/tcllib/
RUN cd /tools/varank/current/tcllib; \
	tar -xvf tcllib.tar.gz
RUN cd /tools/varank/current/ExtAnn; \
	rm -f *; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/ACMG.tsv; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/ClinGenAnnotations.tsv; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/DDG2P.tsv; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/ExAC.CNV-Zscore.annotations.tsv; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/GeneIntolerance-Zscore.annotations.tsv; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/gnomAD.LOEUF.pLI.annotations.tsv; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/HI.tsv; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/morbidGenes.tsv; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/morbidGenesCandidates.tsv; \
	ln -s /STARK/config/variantannotation/varank/configfiles/extanns/default/results.txt; \
	ln -s /STARK/databases/OMIMannotations/current/OMIM-1-annotations.tsv.gz /tools/varank/current/ExtAnn/OMIM-1-annotations.tsv.gz; \
	ln -s /STARK/databases/OMIMannotations/current/OMIM-2-annotations.tsv.gz /tools/varank/current/ExtAnn/OMIM-2-annotations.tsv.gz
RUN cd /tools/varank/current/Annotations_Exomiser; \
	ln -s /STARK/databases/Annotations_Exomiser/current/ 2007

WORKDIR /tools/varank/current/sources/
RUN sed -i -e 's/1.4.3/1.5.2/g' VaRank-config.tcl
RUN sed -i -e '/global L_hpo/a set auto_path [linsert $::auto_path 0 /tools/varank/current/tcllib]' VaRank-main.tcl

ENV VARANK=/tools/varank/current/

#bcftools
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
RUN tar -xvf bcftools-1.19.tar.bz2
WORKDIR /tmp/bcftools-1.19
RUN ./configure --prefix=/tools/
RUN make
RUN make install

# Alamut-Batch
RUN mkdir /tools/alamut_batch
WORKDIR /tools/alamut_batch
RUN wget https://downloads.interactive-biosoftware.com/download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11
RUN tar -xvf download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11; \
	rm -rf download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11; \
	ln -s alamut-batch-standalone-1.11 current
WORKDIR /tools/alamut_batch/current
RUN mv ancillary ancillary.old

COPY src/ancillary.tar.gz /tools/alamut_batch/current/ancillary/
WORKDIR /tools/alamut_batch/current/ancillary
RUN tar -xvf ancillary.tar.gz
RUN rm -rf ancillary.tar.gz

WORKDIR /tools/alamut_batch/current
RUN chmod -R 777 ancillary
RUN rm -rf alamut-batch.ini
RUN ln -s /STARK/config/variantannotation/varank/alamut-batch-license/alamut-batch.ini /tools/alamut_batch/current/alamut-batch.ini

ENV ALAMUT=/tools/alamut_batch/current/
ENV PYTHONDONTWRITEBYTECODE 1

COPY bin/*.py /tools/varank/current/scripts/
RUN chmod -R 777 /tools/

RUN echo "source activate common" >> ~/.bashrc
RUN echo "alias ll='ls -lra'" >> ~/.bashrc
RUN echo "alias varank='python3.12 /tools/varank/current/scripts/__main__.py'" >> ~/.bashrc


WORKDIR /