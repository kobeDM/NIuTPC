// -*- C++ -*-
/*!
 * @file 
 * @brief
 * @date
 * @author
 *
 */

#ifndef NIUTPCMONITOR_H
#define NIUTPCMONITOR_H

#include "DaqComponentBase.h"

#include <arpa/inet.h> // for ntohl()

#include "dc_event.h" //
#include <vector>

////////// ROOT Include files //////////
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TText.h"
#include "TAttText.h"
#include "TNamed.h"

#include "TString.h"
#include "math.h"

#include "TH2.h"
#include "TAxis.h"
#include "TVirtualFFT.h"

#include "TF1.h"

static const int N_BOARD = 4;

static const int HEADER_SIZE = 4*4;
static const int DATA_SIZE = 1024*1024;
static const int N_CHANNEL = 64;
static const int N_STRIP = 32;
static const int N_SAMPLE = 4000;

//2.5MHz -->400nsec/bin 
static const double SAMPLE_RATE = 2500.0;  //kHz

//FFT LowPassFilter  80kHz
static const double LowPassFilter = 80.0;

#include "SampleData.h"

using namespace RTC;

class NIuTPCMonitor
    : public DAQMW::DaqComponentBase
{
public:
    NIuTPCMonitor(RTC::Manager* manager);
    ~NIuTPCMonitor();

    // The initialize action (on CREATED->ALIVE transition)
    // former rtc_init_entry()
    virtual RTC::ReturnCode_t onInitialize();

    // The execution action that is invoked periodically
    // former rtc_active_do()
    virtual RTC::ReturnCode_t onExecute(RTC::UniqueId ec_id);

private:
    TimedOctetSeq          m_in_data;
    InPort<TimedOctetSeq>  m_InPort;

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
    int reset_InPort();

    unsigned int read_InPort();
    //int online_analyze();
    int decode_data(const unsigned char* mydata);
    int fill_data(const unsigned char* mydata, const int size);
    void do_analysis();

    BufferStatus m_in_status;

    int      m_monitor_update_rate;
    const static unsigned int DATA_BUF_SIZE = 1024*1024;
    unsigned char m_recv_data[DATA_BUF_SIZE];
    unsigned int  m_event_byte_size;

    std::vector<dc_event> m_dcevt_coll;
    int       m_ndc;

    ////////////////////////////////////////////////
    //////////////   user   ////////////////////////
    ////////////////////////////////////////////////	

    bool m_debug;
    bool m_set_noise;

		//////////////ROOT declaration//////////////////
		int m_adc_ch;
		int m_adc_xmin;
		int m_adc_xmax;
		double m_mv_ymin;
		double m_mv_ymax;
		double m_pedestal_window;

		int m_prev_evtid;

		// CANVAS
		TCanvas* m_c[N_BOARD];
		TCanvas* m_c_noise;
		//histo
		TH1D* m_h_adc[N_BOARD][N_CHANNEL];
		TH1D* m_h_adc_ped[N_BOARD][N_CHANNEL];
		TH2D* m_h_strips[N_BOARD];

		TH1D* m_h_pedestal[N_BOARD];
		TH1D* m_h_noise[N_BOARD];

		//histo
		TH1D* m_h_adc_HG[N_BOARD][N_STRIP];
		TH1D* m_h_adc_LG[N_BOARD][N_STRIP];
		TH1D* m_h_adc_ped_HG[N_BOARD][N_STRIP];
		TH1D* m_h_adc_ped_LG[N_BOARD][N_STRIP];
		TH2D* m_h_strips_HG[N_BOARD];
		TH2D* m_h_strips_LG[N_BOARD];

		TH1D* m_h_pedestal_HG[N_BOARD];
		TH1D* m_h_pedestal_LG[N_BOARD];
		TH1D* m_h_noise_HG[N_BOARD];
		TH1D* m_h_noise_LG[N_BOARD];

};


extern "C"
{
	void NIuTPCMonitorInit(RTC::Manager* manager);
};

#endif // NIUTPCMONITOR_H
