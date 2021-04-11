#!/usr/bin/python

import os,sys,time

if len(sys.argv) != 3:
    print 'Usage : ',sys.argv[0],' <ip> <size>'
    sys.exit()



ip=str(sys.argv[1])
size=int(sys.argv[2])

size_h = str(format(size,'04x'))
size_h0=size_h[0:2]
size_h1=size_h[2:]

print ip,size,size_h,size_h0,size_h1

os.system("echo \"wrb 0x06 0x%s\n quit\" | python rbcpshell.py %s 4660" % (size_h0,ip))
time.sleep(1)
os.system("echo \"wrb 0x07 0x%s\n quit\" | python rbcpshell.py %s 4660" % (size_h1,ip))
time.sleep(1)
os.system("echo \"rd 0x06 2\n quit\" | python rbcpshell.py %s 4660" % (ip))



