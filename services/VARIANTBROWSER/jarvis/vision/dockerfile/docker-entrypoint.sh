#!/bin/bash

# configure hostname
SERVER_NAME_FILE=/etc/httpd/conf.d/servername.conf
if [ ! -f "$SERVER_NAME_FILE" ]; then
  echo "ServerName $(hostname)" > $SERVER_NAME_FILE
fi


# clean up if container is being restarted
rm -rf /var/run/httpd/*


# run apache
exec httpd -DFOREGROUND
