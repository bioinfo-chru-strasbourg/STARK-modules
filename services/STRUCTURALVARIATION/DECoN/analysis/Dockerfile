
##############################################################
# Dockerfile Version:   0.9
# Software:             DECoN
# Software Version:     1.0.2
# Software Website:     https://github.com/RahmanTeam/DECoN
# Licence:              MIT License
# Description:          DECoN (Detection of Exon Copy Number) is a sensitive and specific tool for detection of exon copy number variants (exon CNVs) in targeted sequencing data
# Usage:                docker run --rm decon:version <RUN> OUTPUT>
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
    Software="DECoN" \
	Version="0.9" \
	Website="https://github.com/RahmanTeam/DECoN" \
	Description="DECoN (Detection of Exon Copy Number) is a sensitive and specific tool for detection of exon copy number variants (exon CNVs) in targeted sequencing data" \
	License="MIT License" \
	Usage="docker run --rm [-v host:container] decon:version <RUN> OUTPUT>"



##############
# PARAMETERS #
##############

ENV TOOLS=/STARK/tools
ENV DATA=/STARK/data
ENV TOOL=/tool
ENV DATABASES=/STARK/databases
ENV YUM_INSTALL="wget gcc-gfortran zlib-devel"
ENV YUM_REMOVE="wget "




###############
# YUM INSTALL #
###############

#RUN yum install -y $YUM_INSTALL ;

# add webtatic for php55
#COPY conf/webtatic.repo /etc/yum.repos.d/webtatic.repo

# install necessary packages
RUN yum upgrade -y && \
	yum updateinfo -y && \
	yum install -y $YUM_INSTALL && \
    yum clean all



################
# DEPENDENCIES #
################




#####
# R #
#####

ENV TOOL_NAME=R
ENV TOOL_VERSION=3.1.2
ENV TARBALL_LOCATION=https://cran.r-project.org/src/base/R-3
ENV TARBALL=R-$TOOL_VERSION.tar.gz
ENV TARBALL_FOLDER=R-$TOOL_VERSION
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL (add webtatic for php55)
RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/Rinstall && \
    cd $TOOLS/$TOOL_NAME/$TOOL_VERSION/Rinstall && \
    wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cd $TARBALL_FOLDER && \
    ./configure --with-readline=no --with-x=no --prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION && \
    make -j $THREADS && \
    make -j $THREADS install && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    wget https://cran.r-project.org/src/contrib/Rcpp_1.0.5.tar.gz && \
    R CMD INSTALL Rcpp_1.0.5.tar.gz && \
    cd $TOOLS/$TOOL_NAME/$TOOL_VERSION && \
    rm -rf Rinstall





#########
# DECoN #
#########


ENV TOOL_NAME=DECoN
ENV TOOL_VERSION=1.0.2
ENV TARBALL_LOCATION=https://github.com/RahmanTeam/DECoN/archive
ENV TARBALL=v$TOOL_VERSION.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME-$TOOL_VERSION
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL (add webtatic for php55)
RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/DECoNinstall && \
    cd $TOOLS/$TOOL_NAME/$TOOL_VERSION/DECoNinstall && \
    wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cp -prf DECoN-1.0.2/Linux/* $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/packrat/src/VGAM $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/packrat/src/scales $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/packrat/src/ggplot2 && \
    wget https://cran.r-project.org/src/contrib/Archive/VGAM/VGAM_0.9-8.tar.gz -O $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/packrat/src/VGAM/VGAM_0.9-8.tar.gz && \
    wget https://cran.r-project.org/src/contrib/Archive/scales/scales_0.2.4.tar.gz -O $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/packrat/src/scales/scales_0.2.4.tar.gz && \
    wget https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_1.0.1.tar.gz -O $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/packrat/src/ggplot2/ggplot2_1.0.1.tar.gz && \
    cd $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
    if [ -w ".Rprofile" ]; then rm .Rprofile; fi && \
    $TOOLS/R/3.1.2/bin/Rscript sessionInfo.R --bootstrap-packrat > setup.log 2>&1 && \
    cp packrat/packrat_source/.Rprofile ./ && \
    ln -s $TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    ln -fs $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOL && \
    rm -rf $TOOLS/$TOOL_NAME/$TOOL_VERSION/DECoNinstall $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/packrat/src



ADD ./scripts/* $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/



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
ENTRYPOINT [ "/tool/bin/DECoN.launch.sh" ]
