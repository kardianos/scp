#!/bin/bash
cd /home/d/code/scp/v102/perm/cas
until [ -f out/x12.done ]; do sleep 30; done
./run_x3c.sh
