// -*- C++ -*-
/*!
 * @file 
 * @brief
 * @date
 * @author
 *
 */

#ifndef NIUTPCREADER0_H
#define NIUTPCREADER0_H

#include <sys/epoll.h>
#include <daqmw/Sock.h>
#include <err.h>

#include "DaqComponentBase.h"
#include "logUtil.h"


static const int N_BOARD=1;

using namespace RTC;

struct module_info {
    DAQMW::Sock Sock;
    std::string ip_address;
    int         port;
    int         module_num;
    unsigned char buf[1024*1024]; // 1024kB buffer
};

class NIuTPCReader0
    : public DAQMW::DaqComponentBase
{
public:
    NIuTPCReader0(RTC::Manager* manager);
    ~NIuTPCReader0();

    // The initialize action (on CREATED->ALIVE transition)
    // former rtc_init_entry()
    virtual RTC::ReturnCode_t onInitialize();

    // The execution action that is invoked periodically
    // former rtc_active_do()
    virtual RTC::ReturnCode_t onExecute(RTC::UniqueId ec_id);

private:
    TimedOctetSeq          m_out_data;
    OutPort<TimedOctetSeq> m_OutPort;

private:
    int daq_dummy();
    int daq_configure();
    int daq_unconfigure();
    int daq_start();
    int daq_run();
    int daq_stop();
    int daq_pause();
    int daq_resume();

    int parse_params(::NVList* list);
    int set_data(unsigned char *buf, unsigned int data_byte_size);
    int write_OutPort();
    int set_window_size(std::string ip_address, int window_size);
    int set_packet_id(std::string ip_address, int module_num);
    int set_delay(std::string ip_address, int delay);
    int set_thinning(std::string ip_address, int thinning);    

    DAQMW::Sock* m_sock;               /// socket for data server

    unsigned int  m_recv_byte_size;
    static const int LAR_TPC_HEADER_BYTE_SIZE    = 4*4;
    static const int LAR_TPC_ONE_EVENT_BYTE_SIZE = 2;
    static const int LAR_TPC_N_CHANNEL           = 64;
    unsigned int  m_window_size;
    unsigned int  m_delay;
    unsigned int  m_thinning;  
    unsigned int  m_read_byte_size;
    int m_epfd;
    struct epoll_event *m_ev_ret;
    bool          m_set_registers;

    BufferStatus m_out_status;

    int m_srcPort;                        /// Port No. of data server
    std::string m_srcAddr;                /// IP addr. of data server
    
    std::vector <struct module_info> m_module_list;

    bool m_debug;
};


extern "C"
{
    void NIuTPCReader0Init(RTC::Manager* manager);
};

#endif // NIUTPCREADER0_H
