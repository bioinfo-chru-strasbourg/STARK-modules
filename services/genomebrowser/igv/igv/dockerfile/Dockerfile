#FROM ubuntu:latest
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive

ARG IGV_WEBAPP_RELEASE=latest

ADD https://raw.githubusercontent.com/dceoy/print-github-tags/master/print-github-tags /usr/local/bin/print-github-tags

RUN set -e \
      && ln -sf bash /bin/sh

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates curl git nginx npm \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && npm install -g yarn

RUN set -eo pipefail \
      && chmod +x /usr/local/bin/print-github-tags \
      && if [ "$IGV_WEBAPP_RELEASE" == "latest" ]; then IGV_WEBAPP_URL=$(print-github-tags --latest --tar igvteam/igv-webapp); else IGV_WEBAPP_URL="https://github.com/igvteam/igv-webapp/archive/v"$IGV_WEBAPP_RELEASE".tar.gz"; fi \
      && curl -SL $IGV_WEBAPP_URL -o /tmp/igv-webapp.tar.gz \
      && tar xvf /tmp/igv-webapp.tar.gz -C /opt --remove-files \
      && mv /opt/igv-webapp-* /opt/igv-webapp \
      && cd /opt/igv-webapp \
      && yarn install \
      && yarn run build

RUN set -e \
      && rm -rf /var/www/html \
      && ln -s /opt/igv-webapp/dist /var/www/html \
      && ln -sf /dev/stdout /var/log/nginx/access.log \
      && ln -sf /dev/stderr /var/log/nginx/error.log \
      && sed -i 's/listen \[::\]:80/#listen [::]:80/g' /etc/nginx/sites-enabled/default

EXPOSE 80

ENTRYPOINT ["/usr/sbin/nginx"]
CMD ["-g", "daemon off;"]