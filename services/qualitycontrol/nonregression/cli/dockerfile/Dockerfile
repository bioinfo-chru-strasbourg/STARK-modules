FROM stark/stark:0.9.18.3

###############################
### Container of nonregression
###############################

RUN yum update -y \
	&& yum install -y \
	figlet \
	git \
	ruby


#set up env 
### FIGLET
WORKDIR /
RUN git clone https://github.com/xero/figlet-fonts.git

### LOLCAT
RUN wget https://github.com/busyloop/lolcat/archive/master.zip
RUN unzip master.zip && cd /lolcat-master/bin
RUN gem install lolcat

## nonregressionapp
COPY . /

CMD ["/bin/bash"]