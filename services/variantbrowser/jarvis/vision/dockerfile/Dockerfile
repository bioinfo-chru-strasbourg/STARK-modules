
##############################################################
# Dockerfile Version:   1.0.2
# Software:             VISION
# Software Version:     1.0.2b
# Software Website:     https://github.com/bioinfo-chru-strasbourg/VISION
# Licence:              GNU Affero General Public License (AGPL)
# Description:          VISION
# Usage:                docker run --rm -d -p HOST_PORT:80 [--name vision -v HOST_CONFIG_FOLDER:/var/www/html/config] vision:version
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
# FROM #
########

FROM centos:7
MAINTAINER Antony Le Bechec <antony.lebechec@gmail.com>
LABEL Software="VISION" \
	Version="1.0.3" \
	Website="https://github.com/bioinfo-chru-strasbourg/VISION/" \
	Description="VISION" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run --rm -d -p HOST_PORT:80 [-v HOST_CONFIG_FOLDER:/var/www/html/config] vision:version"


##############
# PARAMETERS #
##############

ENV TOOLS=/home/TOOLS/tools
ENV DATA=/data
ENV TOOL=/tool
ENV DATABASES=/databases
ENV YUM_INSTALL="gcc bc make wget perl-Switch perl-Digest-MD5 perl-Data-Dumper which zlib-devel zlib 	zlib2-devel zlib2 	bzip2-devel bzip2 	lzma-devel lzma 	xz-devel xz 	ncurses-devel 	unzip"
ENV YUM_REMOVE="zlib-devel zlib2-devel bzip2-devel lzma-devel xz-devel ncurses-devel unzip gcc"




##########
# JARVIS #
##########

ENV TOOL_NAME=VISION
ENV TOOL_VERSION="1.0.3"
ENV TOOL_TARBALL="$TOOL_VERSION.tar.gz"
ENV TARBALL_LOCATION=https://github.com/bioinfo-chru-strasbourg/VISION/archive/refs/heads/
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL (add webtatic for php55)
RUN yum install -y wget && \
	wget --no-check-certificate $TARBALL_LOCATION/$TOOL_TARBALL && \
    tar xf $TOOL_TARBALL && \
    rm -rf $TOOL_TARBALL && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/ && \
    cp $TOOL_NAME-$TOOL_VERSION/* $TOOLS/$TOOL_NAME/$TOOL_VERSION/ -R && \
    rm -rf $TOOL_NAME-$TOOL_VERSION && \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    chmod 0775 $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current -R && \
	mkdir -p $DATABASES && \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ /tool && \
	cp $TOOLS/$TOOL_NAME/$TOOL_VERSION/conf/webtatic.repo /etc/yum.repos.d/webtatic.repo && \
	chmod -R 0777 $TOOLS



###############
# YUM INSTALL #
###############

#RUN yum install -y $YUM_INSTALL ;

# add webtatic for php55
#COPY conf/webtatic.repo /etc/yum.repos.d/webtatic.repo

# install necessary packages
RUN yum upgrade -y && \
	yum updateinfo -y && \
	yum install -y \
	autoconf \
	curl \
	curl-devel \
	cyrus-sasl-lib \
	gcc \
	g++ \
	httpd \
	libc-dev \
	libedit2 \
	libxml2 \
	make \
	php55w \
	php55w-gd \
	php55w-ldap \
	php55w-mysql \
	php55w-pdo \
	php55w-pecl-apcu \
	php55w-pgsql \
	pkg-config \
	which \
	xz-utils \
	openssl-devel \
	psmisc \
	postfix \
	ImageMagick \
	ImageMagick-devel \
	bc \
	make \
	wget \
	ghostscript \
	libgomp \
	java \
	enscript \
	perl-Switch \
	perl-Digest-MD5 \
	perl-Data-Dumper \
	which \
	rsync \
	wine \
	zlib-devel zlib \
    zlib2-devel zlib2 \
    bzip2-devel bzip2 \
    lzma-devel lzma \
    xz-devel xz \
    ncurses-devel \
	$YUM_INSTALL \
    && \
    yum clean all



################
# DEPENDENCIES #
################


##########
# HTSLIB #
##########

ENV TOOL_NAME=htslib
ENV TOOL_VERSION=1.12
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cd $TOOL_NAME-$TOOL_VERSION && \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install && \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current && \
    cd ../ && \
    rm -rf $TOOL_NAME-$TOOL_VERSION



############
# BCFTOOLS #
############

ENV TOOL_NAME=bcftools
ENV TOOL_VERSION=1.12
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL && \
    rm -rf $TARBALL && \
    cd $TOOL_NAME-$TOOL_VERSION && \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install && \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current && \
    cd ../ && \
    rm -rf $TOOL_NAME-$TOOL_VERSION



###########
# HOWARD #
###########

ENV DATABASES=/databases
ENV TOOL_NAME=howard
ENV TOOL_VERSION="0.9.15.6"
ENV TOOL_TARBALL="$TOOL_VERSION.tar.gz"
ENV TARBALL_LOCATION=https://github.com/bioinfo-chru-strasbourg/howard/archive/refs/heads/
ENV TOOL_DATABASE_FOLDER=/home/TOOLS/databases
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN wget $TARBALL_LOCATION/$TOOL_TARBALL && \
    tar xf $TOOL_TARBALL && \
    rm -rf $TOOL_TARBALL && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/ && \
    cp $TOOL_NAME-$TOOL_VERSION/* $TOOLS/$TOOL_NAME/$TOOL_VERSION/ -R && \
    rm -rf $TOOL_NAME-$TOOL_VERSION  && \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
    chmod 0775 $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current -R && \
	mkdir -p $DATABASES && \
	ln -s $DATABASES $TOOL_DATABASE_FOLDER && \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current && \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ /tool && \
	ln -s $TOOLS/ /tools



#################
# JARVIS/VISION #
#################

ENV TOOL_NAME=VISION
ENV TOOL_VERSION="1.0.3"
ENV TOOL_TARBALL="$TOOL_VERSION.tar.gz"
ENV TARBALL_LOCATION=https://github.com/bioinfo-chru-strasbourg/VISION/archive/refs/heads/
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
ADD ./html /var/www/html

RUN cp -R $TOOLS/$TOOL_NAME/$TOOL_VERSION/html/* /var/www/html/ && \
	cp $TOOLS/$TOOL_NAME/$TOOL_VERSION/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh && \
	cp $TOOLS/$TOOL_NAME/$TOOL_VERSION/conf/httpd.conf /etc/httpd/conf/httpd.conf && \
	cp $TOOLS/$TOOL_NAME/$TOOL_VERSION/conf/php.ini /etc/php.ini && \
	mkdir -p /var/www/html/tmp && \
	chmod -R 0777 /var/www/html && \
	chmod -R 0777 $TOOLS



######################
# YUM REMOVE & CLEAR #
######################

#RUN yum erase -y $YUM_REMOVE ; yum clean all ;
RUN yum clean all ;



############
# ENDPOINT #
############

EXPOSE 80

WORKDIR /var/www/html
CMD [ "docker-entrypoint.sh" ]