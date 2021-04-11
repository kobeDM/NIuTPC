#!/bin/bash

ip=$1
#ip="192.168.10.16"

cd /home/msgc/LTARS_DAQ/RBCPshell
sleep 1
echo -e "load tp_tc/slow_test.txt \\n quit" | python rbcpshell_v2.py ${ip} 4660
sleep 1
echo -e "wrb 0x1e 1 \n quit" | python rbcpshell_v2.py ${ip} 4660
#sleep 1

