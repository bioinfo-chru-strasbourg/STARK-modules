##########################################################################
# Dockerfile Version:   1.2
# Software:             POOL
# Software Version:     1.2
# Description:          Dockerfile to build POOL analysis container image
##########################################################################

ARG STARK_IMAGE="stark/stark:0.9.18.2"

########
# FROM #
########

FROM $STARK_IMAGE
LABEL maintainer="NICAISE Samuel <samuel.nicaise@chru-strasbourg.fr>" \
    Software="DEJAVU" \
	Version="1.2" \
	Description="POOL"


ENV SRC=/src
RUN yum -y install  epel-release
RUN yum -y install  gcc gcc-c++ \
					bzip2 bzip2-devel \
					zlib2 zlib2-devel \
					xz xz-devel \
					unzip wget make docker
	
#############
#  python3  #
#############
#RUN yum -y install https://centos7.iuscommunity.org/ius-release.rpm
RUN yum -y install https://repo.ius.io/ius-release-el7.rpm
RUN yum -y install python36u python36u-libs python36u-devel python36u-pip
RUN python -m pip install pandas


#############
# Snakemake #
#############
RUN cd /tmp  && \
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  && \
	bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/lib/miniconda3/ && \
	echo 'export PATH="$PATH:/usr/local/lib/miniconda3/bin"' > /etc/profile.d/miniconda.sh && \
	chmod 755 /etc/profile.d/miniconda.sh && \
	export PATH="$PATH:/usr/local/lib/miniconda3/bin" && \
	/usr/local/lib/miniconda3/bin/conda install -y -c bioconda -c conda-forge snakemake && \
	chmod -R 755 /usr/local/lib/miniconda3/ && \
	rm -f /tmp/Miniconda3-latest-Linux-x86_64.sh;

##########
#  java  #
##########
RUN yum -y install java java-1.8.0

##############
#  Bedtools  #
##############
# RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz  && \
	# tar xf bedtools-2.27.1.tar.gz  && \
	# rm -rf bedtools-2.27.1.tar.gz  && \
	# cd bedtools2  && \
	# make prefix=$SRC/bedtools/2.27.1 install  && \
	# cd ../  && \
	# rm -rf bedtools2;
ENV BEDTOOLS=/STARK/tools/bedtools/current/bin


##########
#  GATK  #
##########
# RUN wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef" -O GenomeAnalysisTK-3.8-1-0.tar.bz2 && \
	# tar -xf GenomeAnalysisTK-3.8-1-0.tar.bz2 && \
	# rm -rf GenomeAnalysisTK-3.8-1-0.tar.bz2 && \
	# cd GenomeAnalysisTK-3.8-1-0-gf15c1c3ef && \
	# mkdir -p $SRC/gatk/3.8-1-0/bin && \
	# cp GenomeAnalysisTK.jar $SRC/gatk/3.8-1-0/bin && \
	# cd ../ && \
	# rm -rf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef;


########
# POOL #
########
COPY . /app
#RUN echo -e '{\n\t"bedtools":"/src/bedtools/2.27.1/bin",\n\t"gatk":"/src/gatk/3.8-1-0/bin/GenomeAnalysisTK.jar",\n\t"java":"/usr/bin/java",\n\t"annotsv":"/src/AnnotSV/2.1",\n\t"R":"/usr/bin/R"\n}' > /app/config/config.json

RUN chmod -R +x /app
ENV PATH="/app/bin:$PATH"

WORKDIR /STARK
ENTRYPOINT []

CMD ["/app/bin/pool","--help"]
