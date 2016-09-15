#!/bin/bash

while [ ! -f Trinity.fasta ]; do
        echo $(date)
        free -h | awk 'NR == 2 {print $3}'
        sleep 300
done
echo "Trinity.fasta detected!"
