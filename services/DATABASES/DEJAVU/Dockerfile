
##############################################################
# Dockerfile Version:   0.9
# Software:             DEJAVU
# Software Version:     0.9
# Software Website:     https://github.com/RahmanTeam/DECoN
# Licence:              GNU Affero General Public License (AGPL)
# Description:          DEJAVU
# Usage:                docker run -ti -v DATABASE_FOLDER:/STARK/databases dejavu:version
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

ARG STARK_IMAGE="stark"
ARG STARK_RELEASE="0.9.18.2"
ARG THREADS="1"



########
# FROM #
########

FROM $STARK_IMAGE:$STARK_RELEASE
LABEL maintainer="Antony Le Bechec <antony.lebechec@gmail.com>" \
    Software="DEJAVU" \
	Version="0.9" \
	Website="https://gitlab.bioinfo-diag.fr/Strasbourg/STARK" \
	Description="DEJAVU" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run -ti -v DATABASE_FOLDER:/STARK/databases dejavu:version"



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




#####
# R #
#####



ADD ./scripts/* $TOOL/toolbox/



######################
# YUM REMOVE & CLEAR #
######################

#RUN yum erase -y $YUM_REMOVE ; yum clean all ;
RUN yum clean all ;



############
# ENDPOINT #
############

EXPOSE 80

WORKDIR /STARK/data
ENTRYPOINT [ "/tool/toolbox/dejavu.generate.sh" ]
