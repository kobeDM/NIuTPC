#!/bin/bash

ip=$1
#ip="192.168.13.11"

cd /home/msgc/LTARS_DAQ/RBCPshell
sleep 1
echo -e "load mode/fast_mode.txt \\n quit" | python rbcpshell.py ${ip} 4660
sleep 1
#echo -e "wrb 0x1e 1 \\n quit " | python rbcpshell.py ${ip} 4660
#sleep 1

