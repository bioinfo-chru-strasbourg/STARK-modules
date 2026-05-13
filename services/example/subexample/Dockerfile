FROM almalinux:9
LABEL maintainer="bioinfo@chru-strasbourg.fr" \
	Version="2.0.0" \
	Description="subexample image"

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
	python -m pip install cyvcf2 fastapi "uvicorn[standard]" python-dotenv && \
	python -m pip cache purge && \
	conda clean --all --yes

WORKDIR /app

ARG STARK_MAIN_FOLDER
ENV DOCKER_STARK_MAIN_FOLDER=${STARK_MAIN_FOLDER}
RUN echo "DOCKER_STARK_MAIN_FOLDER=${DOCKER_STARK_MAIN_FOLDER}" > /app/.env
