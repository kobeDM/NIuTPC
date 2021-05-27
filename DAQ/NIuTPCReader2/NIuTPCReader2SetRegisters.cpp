#include <arpa/inet.h>
#include <err.h>
#include <cstdio>
#include <iostream>

#include <daqmw/Sock.h>
#include "NIuTPCReader2.h"
#include "SitcpRbcp.h"

int NIuTPCReader2::set_window_size(std::string ip_address, int window_size)
{
    SitcpRbcp bcp;

    int length  = 2;
    int address = 0x6;
    unsigned char data_buf[length];
    //    sprintf(data_buf,"%x",window_size);
    // below sentence is not correct?
    // Register value for Window size
    // Register value:  0 means 1
    // Register value:  1 means 2
    // Register value: 29 means 30
    //        data_buf[0] = window_size - 1;
    data_buf[1] = window_size & 0xff;
    window_size = window_size >> 8;
    data_buf[0] = window_size & 0xff;

    if (bcp.write_registers(ip_address, address, length, data_buf) < 0) {
        warnx("registers for window_size setting fail");
        return -1;
    }

    return 0;
}

int NIuTPCReader2::set_packet_id(std::string ip_address, int module_num)
{
    SitcpRbcp bcp;

    int length  = 3;
    int address = 0x0a;
    unsigned char data_buf[length];

    data_buf[0] = 0x00;
    data_buf[1] = 0x00;
    data_buf[2] = module_num;

    if (bcp.write_registers(ip_address, address, length, data_buf) < 0) {
        warnx("registers for packet_id setting fail");
        return -1;
    }

    return 0;
}

int NIuTPCReader2::set_delay(std::string ip_address, int delay)
{
    SitcpRbcp bcp;

    int length  = 2;
    int address = 0x08;
    unsigned char data_buf[length];
    
    data_buf[1] = delay & 0xff;
    delay = delay >> 8;
    data_buf[0] = delay & 0xff;


    if (bcp.write_registers(ip_address, address, length, data_buf) < 0) {
        warnx("registers for delay setting fail");
        return -1;
    }

    return 0;
}
int NIuTPCReader2::set_thinning(std::string ip_address, int thinning)
{
    SitcpRbcp bcp;

    int length  = 1;
    int address = 0x0d;
    unsigned char data_buf[length];
    
    data_buf[0] = thinning & 0xff;


    if (bcp.write_registers(ip_address, address, length, data_buf) < 0) {
        warnx("registers for thinning setting fail");
        return -1;
    }

    return 0;
}
