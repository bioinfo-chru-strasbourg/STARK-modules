#!/bin/bash
#only configure once
#[ -f /var/run/slurmctld.startup ] && exit 0

sed -e '/^hosts:/d' -i /etc/nsswitch.conf
echo 'hosts:      files myhostname' >> /etc/nsswitch.conf

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




# wait for cluster
while true
do
	sacctmgr show cluster &>/dev/null
	[ $? -eq 0 ] && break
	sleep 5
done


CLUSTERNAME="$(awk -F= '/^ClusterName=/ {print $2}' /etc/slurm/slurm.conf)"
[ -z "$CLUSTERNAME" ] && echo 'no cluster name' && exit 1

# Add cluster
sacctmgr -vi add cluster "$CLUSTERNAME"
# Add account
sacctmgr -vi add account admin Cluster="$CLUSTERNAME" Description="none" Organization="none"
# Add user
sacctmgr -vi add user root Account=admin DefaultAccount=admin
sacctmgr -vi add user slurm Account=admin DefaultAccount=admin
# Add QOSs:
sacctmgr -i add qos normal
sacctmgr -i add qos high
sacctmgr -i add qos medium
sacctmgr -i add qos low
# Modify user QOS
sacctmgr -i modify account admin set qos=high
sacctmgr -i modify user root,slurm set qos=high

# Load users
if [ ! -s /config/users.sh ]; then
	touch /config/users.sh
fi;
chmod u+x /config/users.sh
/config/users.sh

# Load QOS
if [ ! -s /config/qos.sh ]; then
	touch /config/qos.sh
fi;
chmod u+x /config/qos.sh
/config/qos.sh

# Load cluster.cfg
if [ ! -s /config/cluster.cfg ]; then
	touch /config/cluster.cfg
fi;
sacctmgr -vi load /config/cluster.cfg

if [ "$(hostname -s)" = "mgmtnode" ]
then
	if [ ! -s /etc/slurm/nodes.conf ]
	then
		props="$(slurmd -C | head -1 | sed 's#NodeName=mgmtnode ##g')"
		echo "NodeName=DEFAULT $props Gres=gpu:gtx:3" >> /etc/slurm/nodes.conf

		cat /etc/nodelist | while read name ip4 ip6
		do
			[ ! -z "$ip6" ] && addr="$ip6" || addr="$ip4"
			echo "NodeName=$name NodeAddr=$addr" >> /etc/slurm/nodes.conf
		done
	fi
else
	while [ ! -s /etc/slurm/nodes.conf ]
	do
		sleep 0.25
	done
fi

touch /var/run/slurmctld.startup

exit 0
