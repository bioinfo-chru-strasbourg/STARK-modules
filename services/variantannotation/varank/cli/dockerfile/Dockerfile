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
# varank
RUN yum -y install git
RUN yum -y install glibc-2.17-326.el7_9.i686
RUN yum clean all

#conda
RUN cd /tmp
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/lib/miniconda3/
RUN echo 'export PATH="$PATH:/usr/local/lib/miniconda3/bin"' > /etc/profile.d/miniconda.sh
RUN chmod 755 /etc/profile.d/miniconda.sh
ENV PATH="$PATH:/usr/local/lib/miniconda3/bin"
RUN chmod -R 755 /usr/local/lib/miniconda3/
RUN rm -f /tmp/Miniconda3-latest-Linux-x86_64.sh

#work env
RUN conda install -c intel tcl
RUN conda create -n common python=3.9.7 

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
RUN cd current/ExtAnn; \
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
	ln -s /STARK/databases/OMIMannotations/current/OMIM-1-annotations.tsv.gz /tools/varank/current/ExtAnn/OMIM-1-annotations.tsv.gz; \
	ln -s /STARK/databases/OMIMannotations/current/OMIM-2-annotations.tsv.gz /tools/varank/current/ExtAnn/OMIM-2-annotations.tsv.gz

ENV VARANK=/tools/varank/current/

# Alamut-Batch
RUN mkdir /tools/alamut_batch
WORKDIR /tools/alamut_batch
RUN wget https://downloads.interactive-biosoftware.com/download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11
RUN tar -xvf download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11; \
	rm -rf download.php?product=Alamut-batch-standalone\&platform=\&OS=Linux\&version=1.11; \
	ln -s alamut-batch-standalone-1.11 current
WORKDIR /tools/alamut_batch/current
RUN rm -rf ancillary

COPY src/ancillary.tar.gz /tools/alamut_batch/current/ancillary/
WORKDIR /tools/alamut_batch/current/ancillary
RUN tar -xvf ancillary.tar.gz
RUN rm -rf ancillary.tar.gz
WORKDIR /tools/alamut_batch/current
RUN chmod -R 777 ancillary
RUN rm -rf alamut-batch.ini
RUN ln -s /STARK/config/variantannotation/varank/alamut-batch-license/alamut-batch.ini /tools/alamut_batch/current/alamut-batch.ini
ENV ALAMUT=/tools/alamut_batch/current/

COPY bin/*.py /tools/varank/current/scripts/
RUN chmod -R 777 /tools/

RUN echo "source activate common" > ~/.bashrc
RUN echo "alias ll='ls -l'" > ~/.bashrc
RUN echo "alias varank='python3.9 /tools/varank/current/scripts/__main__.py'" > ~/.bashrc


WORKDIR /