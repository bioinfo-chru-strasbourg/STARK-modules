#!/bin/bash
mkdir -p /service
/usr/local/bin/sacct --allusers -l --json > /service/$(hostname).sacct.json
