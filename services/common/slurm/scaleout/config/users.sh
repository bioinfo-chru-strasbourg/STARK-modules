#!/bin/bash
#

# Create groups:
groupadd -g 1005 admin
groupadd -g 1006 bioinfo
groupadd -g 1007 biology

# Create users
for i in stark fred wilma; do \
    /usr/sbin/useradd -c $i -m -g users $i && \
    echo 'password' | passwd $i --stdin && \
    mkdir -p /home/$i/.ssh && \
    mkdir -p /srv/containers/$i/ && \
    chmod 0700 /srv/containers/$i/ && \
    ssh-keygen -f /home/$i/.ssh/id_ecdsa -N '' -t ecdsa && \
    ssh-keygen -y -f /home/$i/.ssh/id_ecdsa > /home/$i/.ssh/id_ecdsa.pub && \
    cat /home/$i/.ssh/id_ecdsa.pub >> /home/$i/.ssh/authorized_keys && \
    chown -R $i:users /home/$i/.ssh && chmod -R 0700 /home/$i/.ssh; \
    mkdir -p /home/$i/.enroot/{data,cache}; \
done;

# Assign group bionfo to users
for i in stark wilma; do \
    /usr/sbin/usermod -a -G bioinfo stark && \
    chown $i:bioinfo /srv/containers/$i/ && \
    /usr/sbin/usermod -G bioinfo $i; \
done;

# Assign group biology to users
for i in stark fred; do \
    /usr/sbin/usermod -a -G bioinfo stark && \
    chown $i:bioinfo /srv/containers/$i/ && \
    /usr/sbin/usermod -G bioinfo $i; \
done;
