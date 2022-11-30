#!/bin/bash
#

# Create the QOSs:
sacctmgr -i add qos normal
sacctmgr -i add qos high
sacctmgr -i add qos medium
sacctmgr -i add qos low
sacctmgr -i add qos urgent

# QOS priority
sacctmgr -i modify qos normal set Priority=50
sacctmgr -i modify qos high set Priority=100
sacctmgr -i modify qos medium set Priority=50
sacctmgr -i modify qos low set Priority=10
sacctmgr -i modify qos urgent set Priority=1000

# Assign all users, cluster and accounts to the low QOS, and make the low
# use cluster.cfg