#FROM alpine
FROM php:8.0.5-fpm-alpine
RUN mkdir -p /run/nginx/ /var/www/html/validate /var/www/html/auth /var/www/html/dashboard
RUN mkdir -p /var/www/html/slurm
RUN apk update && apk add nginx php8-fpm bash php8-session nginx-debug
RUN unlink /etc/nginx/http.d/default.conf
COPY www.conf /etc/php8/php-fpm.d/www.conf
COPY nginx.conf /etc/nginx/http.d/
COPY validate.php /var/www/html/validate/index.php
COPY auth.php /var/www/html/auth/index.php
COPY dashboard /var/www/html/dashboard
RUN chmod -R 0755 /var/www/html/validate/index.php /var/www/html/auth/index.php /var/www/html/dashboard/
RUN echo "daemon off;" >> /etc/nginx/nginx.conf
RUN echo "error_log syslog:server=unix:/dev/log,nohostname,tag=nginx debug;" >> /etc/nginx/nginx.conf
CMD bash -c 'nginx & php-fpm8 -F'

EXPOSE 8080/tcp
