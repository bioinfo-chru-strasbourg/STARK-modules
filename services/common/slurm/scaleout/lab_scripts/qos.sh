#!/bin/bash
#
# Script to setup for the QOS-based preemption lab

echo "PriorityWeightQOS=1000000" >> /etc/slurm/slurm.conf
/lab_scripts/restart.sh
sleep 3
#Create the QOSs:
sacctmgr -i add qos high
sacctmgr -i add qos medium
sacctmgr -i add qos low
# Assign all users, cluster and accounts to the low QOS, and make the low
# QOS the default for all users, cluster and accounts:
sacctmgr -i modify account admin set qos=high
sacctmgr -i modify account biology set qos=low
sacctmgr -i modify account bioinfo set qos=medium
sacctmgr -i modify user root,fred,slurm,wilma set qos=low
sacctmgr -i modify cluster cluster set qos=low
sacctmgr -i modify account admin set defaultqos=high
sacctmgr -i modify account biology set defaultqos=low
sacctmgr -i modify account bioinfo set defaultqos=medium
sacctmgr -i modify cluster cluster set defaultqos=low

#Assign users to be able to use the QOSs:
sacctmgr -i modify user fred,wilma set qos=+high,+medium,+low
sacctmgr -i modify user fred,wilma set qos=+medium,+low
sacctmgr -i modify user root set qos=+high,+medium,+low

