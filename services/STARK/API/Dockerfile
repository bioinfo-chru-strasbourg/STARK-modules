FROM alpine:3.1

# Update
RUN apk add --update python py-pip docker git make gcc g++

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

# Clean packages
RUN apk del git make gcc g++

EXPOSE  8000
CMD ["python", "/src/API.py", "-p 8000"]
