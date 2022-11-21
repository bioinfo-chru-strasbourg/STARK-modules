#!/bin/bash

# Add hosts in the not crazy slow manner
cat /etc/hosts.nodes >> /etc/hosts

#ensure the systemd cgroup directory exists for enroot
mkdir -p $(awk -F: '$2 ~ /systemd/ {printf "/sys/fs/cgroup/systemd/%s", $3}' /proc/self/cgroup)

#systemd user@.service handles on normal nodes
# for i in stark fred wilma; do
# 	uid=$(id -u $i)
# 	mkdir -m 0700 -p /run/user/$uid
# 	chown $i:users /run/user/$uid
# done
# Load users
if [ ! -s /config/users.sh ]; then
	touch /config/users.sh
fi;
chmod u+x /config/users.sh
/config/users.sh

#start systemd
exec /lib/systemd/systemd --system --log-level=info --crash-reboot --log-target=console
