#!/usr/bin/python

import os,sys

for i in range(64):
    print 'wrb 0x1c 0x%x' % i 
    #print 'wrb 0x1b 0x20'
    print 'wrb 0x1b 0x0'
    print 'wrb 0x1d 0'
    
