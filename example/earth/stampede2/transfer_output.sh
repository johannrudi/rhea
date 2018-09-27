#!/bin/bash

BASEDIR="${1:-.}"
SSH_PORT="8704"
SSH_URL="nordend.ices.utexas.edu"
SSH_DIR="/run/media/johann/RESEARCH/runs/rhea/2018-07_earth"

transfer_file="runs.tgz"
scp -P $SSH_PORT "$transfer_file" ${SSH_URL}:${SSH_DIR}
