#!/bin/bash
#only configure once
#[ -f /var/run/slurmctld.startup ] && exit 0

sed -e '/^hosts:/d' -i /etc/nsswitch.conf
echo 'hosts:      files myhostname' >> /etc/nsswitch.conf

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
