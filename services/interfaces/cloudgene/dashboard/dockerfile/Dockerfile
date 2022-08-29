
##############################################################
# Dockerfile Version:   1.0
# Software:             CloudGene
# Software Version:     2.5.3
# Software Website:     https://github.com/genepi/cloudgene
# Licence:              GNU Affero General Public License (AGPL)
# Description:          CloudGene
# Usage:                docker run -ti -p 8080:80 -v DATA:/STARK/data cloudgene:version
##############################################################

##########
# README #
##########

# Config parameters
#    identify yum packages for installation
#    identify yum packages to remove
#
# Dependecies installation
#    identify tools dependences
#    config each tool
#    write installation procedure for each tools
#
# Tool
#    configure tool
#    write isntallation procedure for the tool
#    add link to current and root tool folder
#
# Workdir / Entrypoint / Cmd
#    configure workdir, endpoint and command
#    /!\ no variables in endpoint



########
# ARGS #
#######

#ARG STARK_IMAGE="stark"
#ARG STARK_RELEASE="0.9.18.2"
ARG STARK_IMAGE="stark/stark:0.9.18.4"
ARG THREADS="1"
ARG CLOUDGENE_RELEASE="2.5.3"


########
# FROM #
########

FROM $STARK_IMAGE
LABEL maintainer="Antony Le Bechec <antony.lebechec@gmail.com>" \
    Software="CloudGene" \
	Version="1.0" \
	Website="https://github.com/genepi/cloudgene" \
	Description="CloudGene" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run -ti -p 8080:80 -v DATA:/STARK/data cloudgene:version"



##############
# PARAMETERS #
##############

ENV TOOLS=/STARK/tools
ENV DATA=/STARK/data
ENV TOOL=/tool
ENV DATABASES=/STARK/databases
ENV YUM_INSTALL=""
ENV YUM_REMOVE=""




###############
# YUM INSTALL #
###############

#RUN yum install -y $YUM_INSTALL ;

# add webtatic for php55
#COPY conf/webtatic.repo /etc/yum.repos.d/webtatic.repo

# install necessary packages
# RUN yum upgrade -y && \
# 	yum updateinfo -y && \
# 	yum install -y $YUM_INSTALL && \
#     yum clean all



################
# DEPENDENCIES #
################


#############
# CloudGene #
#############

# TOOL INFO
ENV TOOL_NAME="cloudgene"
ENV TOOL_VERSION=$CLOUDGENE_RELEASE
ENV TOOL_TARBALL="v$TOOL_VERSION.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/genepi/cloudgene/archive/refs/tags/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# TOOL PARAMETERS
ENV TOOL_JAR=$TOOL_TARBALL
ENV TOOL_JAR_TARGET="VcfStats.jar"

# TOOL INSTALLATION
RUN if [ -z $TOOL_VERSION ]; then TOOL_VERSION_FOLDER=latest; else TOOL_VERSION_FOLDER=$TOOL_VERSION; fi && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin && \
    (cd $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin && \
    curl -s install.cloudgene.io | bash -s $TOOL_VERSION && \
	ln -s $TOOL_VERSION_FOLDER $TOOLS/$TOOL_NAME/current && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/config/apps && \
    ln -s config/apps $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/apps && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/data/tmp && \
    ln -s data/tmp $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/tmp && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/data/workspace && \
    ln -s data/workspace $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/workspace)


############
# VCFSTATS #
############

# # TOOL INFO
# ENV TOOL_NAME="vcfstats"
# ENV TOOL_VERSION="1.2"
# ENV TOOL_TARBALL="VcfStats-assembly-$TOOL_VERSION.jar"
# ENV TOOL_SOURCE_EXTERNAL="https://github.com/biopet/vcfstats/releases/download/v$TOOL_VERSION/$TOOL_TARBALL"
# ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# # TOOL PARAMETERS
# ENV TOOL_JAR=$TOOL_TARBALL
# ENV TOOL_JAR_TARGET="VcfStats.jar"

# # TOOL INSTALLATION
# RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
# 	wget $TOOL_SOURCE_EXTERNAL -O $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/$TOOL_JAR_TARGET && \
# 	ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current;



##########
# Script #
##########

ADD ./scripts/launch.sh $TOOLS/$TOOL_NAME/current/bin/
ADD ./scripts/*.stache $TOOLS/$TOOL_NAME/current/bin/sample/pages/
ADD ./conf/* $TOOLS/$TOOL_NAME/current/bin/config/



######################
# YUM REMOVE & CLEAR #
######################

#RUN yum erase -y $YUM_REMOVE ; yum clean all ;
RUN yum clean all ;



############
# ENDPOINT #
############

# Cloudgene Docker Branding
ENV CLOUDGENE_SERVICE_NAME="Cloudgene Docker"
ENV CLOUDGENE_HELP_PAGE="https://github.com/lukfor/docker-cloudgene"
ENV START_CLOUDGENE="true"
ENV START_HADOOP="true"

EXPOSE 8082

WORKDIR /STARK/data
#ENTRYPOINT [ "/STARK/tools/cloudgene/current/bin/cloudgene", "server" ]
ENTRYPOINT [ "/STARK/tools/cloudgene/current/bin/launch.sh" ]
