#!/bin/bash

ip=$1
id=$2

sleep 1
echo -e "load ./setup/setup1.txt \n quit" | python rbcpshell.py ${ip} 4660
sleep 1
echo -e "load ./setup/setup2.txt \n quit" | python rbcpshell.py ${ip} 4660
sleep 1
echo -e "load ./setup/setup3.txt \n quit" | python rbcpshell.py ${ip} 4660
sleep 1


###echo "load ./setup/setup4.txt\n quit" | python rbcpshell.py ${ip} 4660
##sleep 1
##echo "load ./setup/setup5.txt\n quit" | python rbcpshell.py ${ip} 4660
##sleep 1

echo -e "wrb 0x0a 0x00\n quit" | python rbcpshell.py ${ip} 4660
sleep 1
echo -e "wrb 0x0b 0x00\n quit" | python rbcpshell.py ${ip} 4660
sleep 1
echo -e "wrb 0x0c ${id}\n quit" | python rbcpshell.py ${ip} 4660
sleep 1
echo -e "wrb 0x08 4\n quit" | python rbcpshell.py ${ip} 4660
sleep 1
echo -e "rd 0x0a 3\n quit" | python rbcpshell.py ${ip} 4660


# self trigger (20201204 Satoshi)
sleep 1

# trigger enable : 0 for reset
echo -e "wrb 0x12 0\n quit" | python rbcpshell.py ${ip} 4660
sleep 1

echo -e "wrb 0x12 1\n quit" | python rbcpshell.py ${ip} 4660
sleep 1

# self-trigger threshold for each channe [15:8]
echo -e "wrb 0x13 1\n quit" | python rbcpshell.py ${ip} 4660
sleep 1

# self-trigger threshold for each channe [7:0]
<<<<<<< HEAD
echo -e "wrb 0x14 128\n quit" | python rbcpshell.py ${ip} 4660
# echo -e "wrb 0x14 150\n quit" | python rbcpshell.py ${ip} 4660
# echo -e "wrb 0x14 0\n quit" | python rbcpshell.py ${ip} 4660
=======
# echo -e "wrb 0x14 128\n quit" | python rbcpshell.py ${ip} 4660
echo -e "wrb 0x14 0\n quit" | python rbcpshell.py ${ip} 4660
>>>>>>> e6ca5a6d329073262616eca4059112301fbedca0
sleep 1

# the number of channels used for self-trigger
echo -e "wrb 0x15 3\n quit" | python rbcpshell.py ${ip} 4660

