FROM almalinux:8

# Update
RUN yum update -y && yum install -y epel-release && yum install -y python3 python3-pip docker git make gcc gcc-c++ && yum clean all

# Install app dependencies
RUN pip3 install Flask

# Task Spooler
RUN git clone https://github.com/thomaspreece/task-spooler.git && \
    cd task-spooler && \
    make && \
    make install && \
    mkdir /ts-tmp && \
    rm -rf task-spooler ;

# Bundle app source
COPY API.py /src/API.py

# YUM clean
RUN yum erase -y git make gcc g++ && yum clean all ;

EXPOSE  8000
CMD ["python3", "/src/API.py", "-p 8000"]