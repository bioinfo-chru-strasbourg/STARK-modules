# Dockerfile
FROM quay.io/biocontainers/gemini:0.30.2--py27h6395249_2

RUN mkdir -p /app/config /app/gemini_data && \
    ln -s /app/config/ /usr/local/share/gemini && \
    ln -s /app/gemini_data/ /usr/local/share/gemini/gemini_data

ADD config/* /app/config/

ENTRYPOINT [ "gemini" ]
