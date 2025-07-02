FROM almalinux:8

EXPOSE 5000
WORKDIR /app

RUN yum install python3 python3-pip git -y && \
    git clone https://github.com/igvteam/igv.js-flask.git . && \
    rm -rf /app/test /app/test_IGV.py && \
    rm -rf /app/igvjs/static/data/public/gstt1_sample.bam /app/igvjs/static/data/public/gstt1_sample.bam.bai && \
    rm -rf /app/igvjs/static/data/public/test.bedGraph /app/igvjs/static/data/public/test.bw && \
    pip3 install -r requirements.txt;

COPY help.html /app/igvjs/static/help.html

CMD ["python3", "/app/run.py"]
