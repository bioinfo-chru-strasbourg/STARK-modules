#!/usr/bin/bash

echo "$@" | systemd-cat -t scrun-trace

exec /usr/local/sbin/scrun "$@"
