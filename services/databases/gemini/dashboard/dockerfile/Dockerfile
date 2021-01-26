# Dockerfile
FROM php:7.4-apache

ENV COMPOSER_ALLOW_SUPERUSER=1

EXPOSE 80
WORKDIR /app

RUN apt-get update -y && \
    apt-get install -y git libicu-dev unzip zip python wget gcc && \
    curl -sS https://getcomposer.org/installer | php -- --install-dir=/usr/local/bin --filename=composer && \
	apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /STARK/tools/gemini/latest && \
    cd /STARK/tools/gemini/latest && \
    wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py && \
    python gemini_install.py /STARK/tools/gemini/latest /STARK/tools/gemini/latest/data --nodata && \
    chmod a+x -R /STARK/tools/gemini/latest/*

ENV PATH=$PATH:/STARK/tools/gemini/latest/bin

# PHP Extensions
RUN docker-php-ext-install -j$(nproc) opcache pdo_mysql
ADD conf/php.ini /usr/local/etc/php/conf.d/app.ini

# Apache
#ADD www/errors /errors
RUN ln -s /app/errors/ /errors
RUN a2enmod rewrite remoteip
ADD conf/vhost.conf /etc/apache2/sites-available/000-default.conf
ADD conf/apache.conf /etc/apache2/conf-available/z-app.conf
RUN a2enconf z-app
ADD www /app
RUN chmod 0755 /app/* -R
