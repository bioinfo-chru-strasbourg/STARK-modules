FROM almalinux:9
LABEL maintainer="bioinfo@chru-strasbourg.fr" \
	Version="2.0.0" \
	Description="POOL image"

ENV LC_CTYPE=en_US.UTF-8

RUN dnf -y update && \
	dnf -y install python3.12 wget

RUN dnf -y install dnf-plugins-core && \
	dnf config-manager --add-repo=https://download.docker.com/linux/centos/docker-ce.repo && \
	dnf -y install docker-ce docker-ce-cli containerd.io

WORKDIR /tmp
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
	bash Miniforge3-Linux-x86_64.sh -b -p /opt/conda && \
	rm -rf Miniforge3-Linux-x86_64.sh && \
	chmod -R 777 /opt/conda
ENV PATH="/opt/conda/bin:$PATH"

# htslib-tools contains bgzip and tabix
RUN dnf install -y epel-release && \
	dnf -y install bcftools htslib htslib-tools

COPY . /app
RUN chmod -R 777 /app

WORKDIR /app
RUN /opt/conda/bin/conda init --all && \
	/opt/conda/bin/conda create -n py312 python=3.12 -y && \
	source /opt/conda/etc/profile.d/conda.sh && \
	conda activate py312 && \
	python -m ensurepip --upgrade && \
	python -m pip install --upgrade pip && \
	python -m pip install cyvcf2 && \
	python -m pip cache purge && \
	conda clean --all --yes

# mosdepth
# RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.12/mosdepth && \
# 	chmod +x mosdepth && \
# 	mv mosdepth /usr/local/bin/

# pandepth
# WORKDIR /tmp
# RUN wget -c https://github.com/HuiyangYu/PanDepth/releases/download/v2.26/PanDepth-2.26-Linux-x86_64.tar.gz && \
# 	tar zvxf PanDepth-2.26-Linux-x86_64.tar.gz && \
# 	rm -rf PanDepth-2.26-Linux-x86_64.tar.gz && \
# 	cd PanDepth-2.26-Linux-x86_64 && \
# 	chmod +x pandepth && \
# 	mv pandepth /usr/local/bin/

# bedtools
# RUN dnf -y install gcc make zlib-devel bzip2-devel xz-devel gcc-c++ && \
# 	wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
# 	tar -zxvf bedtools-2.31.1.tar.gz && \
# 	rm -rf bedtools-2.31.1.tar.gz && \
# 	cd bedtools2 &&	\
# 	make && \
# 	mv bin/* /usr/local/bin/

WORKDIR /app

ARG STARK_MAIN_FOLDER
ENV DOCKER_STARK_MAIN_FOLDER=${STARK_MAIN_FOLDER}
RUN echo "DOCKER_STARK_MAIN_FOLDER=${DOCKER_STARK_MAIN_FOLDER}" > /app/.env

# WORKDIR /STARK
# ENTRYPOINT []

# CMD ["/app/bin/pool","--help"]
