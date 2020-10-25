FROM centos:7

# Update
RUN yum update -y && yum install -y epel-release && yum install -y python python-pip docker git make gcc g++ && yum clean all
#RUN yum install -y python-pip
#RUN yum install -y docker


# Install app dependencies
RUN pip install Flask


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
CMD ["python", "/src/API.py", "-p 8000"]
