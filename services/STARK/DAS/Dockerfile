# Dockerfile
FROM centos:7

EXPOSE 5000
WORKDIR /app

RUN yum install python3 python3-pip git -y && \
    git clone https://github.com/igvteam/igv.js-flask.git . && \
    pip3 install -r requirements.txt;

COPY help.html /app/igvjs/static/help.html

CMD ["python3", "/app/run.py"]
