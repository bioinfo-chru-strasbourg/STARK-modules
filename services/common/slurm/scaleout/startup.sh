#!/bin/bash


# Nodes
if [ ! -s /config/nodelist ] || [ ! -s /config/hosts.nodes ] || [ ! -s /config/hosts ]; then
	touch /config/nodelist /config/hosts.nodes /config/hosts
fi;

if [ -s /etc/hosts ]; then
	cat /etc/hosts /config/hosts | sort -u >> /config/hosts
fi;

if [ ! -z "$NODE_NAME" ]; then
	if ! grep -Fxq "$NODE_NAME $NODE_IP4 $NODE_IP6" /config/nodelist; then
		echo -e "$NODE_NAME $NODE_IP4 $NODE_IP6" >> /config/nodelist
	fi
	#echo -e "$NODE_NAME $NODE_IP4 $NODE_IP6" >> /config/nodelist
	if ! grep -Fxq ""$NODE_IP4 $NODE_NAME"" /config/hosts.nodes; then
		echo -e ""$NODE_IP4 $NODE_NAME"" >> /config/hosts.nodes
	fi
	if ! grep -Fxq ""$NODE_IP6 $NODE_NAME"" /config/hosts.nodes; then
		echo -e ""$NODE_IP6 $NODE_NAME"" >> /config/hosts.nodes
	fi
	# echo -e "$NODE_IP4 $NODE_NAME" >> /config/hosts.nodes
	# echo -e "$NODE_IP6 $NODE_NAME" >> /config/hosts.nodes
	if ! grep -Fxq ""$NODE_IP4 $NODE_NAME"" /config/hosts; then
		echo -e ""$NODE_IP4 $NODE_NAME"" >> /config/hosts
	fi
	if ! grep -Fxq ""$NODE_IP6 $NODE_NAME"" /config/hosts; then
		echo -e ""$NODE_IP6 $NODE_NAME"" >> /config/hosts
	fi
fi;

awk -i inplace '!seen[$0]++'  /config/nodelist
awk -i inplace '!seen[$0]++'  /config/hosts.nodes
awk -i inplace '!seen[$0]++'  /config/hosts

ln -sfn /config/nodelist /etc/nodelist
ln -sfn /config/hosts.nodes /etc/hosts.nodes
ln -sfn /config/hosts /etc/hosts

# yes | cp /config/nodelist /etc/nodelist
# yes | cp /config/hosts.nodes /etc/hosts.nodes

# Add hosts in the not crazy slow manner
#cat /config/hosts.nodes >> /etc/hosts



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
