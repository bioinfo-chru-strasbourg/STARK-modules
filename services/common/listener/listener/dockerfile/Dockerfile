FROM centos:7
MAINTAINER Victor Grentzinger <victor.grentzinger@gmail.com>

##################################################################################
# Dockerfile Version:   1.0
# Software:             listener
# Software Version:     1.0
# Description:          Dockerfile to build Microservices Listener container image
##################################################################################


#######
# YUM #
#######
RUN yum update -y
RUN yum -y install python2 epel-release yum-utils device-mapper-persistent-data lvm2

#############
# DOCKER-CE #
#############
RUN yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
RUN yum -y install docker-ce
RUN usermod -aG docker $(whoami)
RUN systemctl enable docker.service

############
# LISTENER #
############

COPY . /app

RUN chmod -R +x /app
ENV PATH="/app:$PATH"

WORKDIR /app

CMD ["python","listener.py","--help"]
