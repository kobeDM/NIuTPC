// -*- C++ -*-
/*!
 * @file
 * @brief
 * @date
 * @author
 *
 */

#include "NIuTPCMonitor.h"
#include <assert.h>
#include <vector>
#include <numeric>

using DAQMW::FatalType::DATAPATH_DISCONNECTED;
using DAQMW::FatalType::INPORT_ERROR;
using DAQMW::FatalType::HEADER_DATA_MISMATCH;
using DAQMW::FatalType::FOOTER_DATA_MISMATCH;
using DAQMW::FatalType::USER_DEFINED_ERROR1;

// Module specification
// Change following items to suit your component's spec.
static const char* niutpcmonitor_spec[] =
    {
        "implementation_id", "NIuTPCMonitor",
        "type_name",         "NIuTPCMonitor",
        "description",       "NIuTPCMonitor component",
        "version",           "1.0",
        "vendor",            "Kazuo Nakayoshi, KEK",
        "category",          "example",
		"activity_type",     "DataFlowComponent",
		"max_instance",      "1",
		"language",          "C++",
		"lang_type",         "compile",
		""
    };

NIuTPCMonitor::NIuTPCMonitor(RTC::Manager* manager)
	: DAQMW::DaqComponentBase(manager),
	m_InPort("niutpcmonitor_in",   m_in_data),
	m_in_status(BUF_SUCCESS),
	m_monitor_update_rate(30),
	m_event_byte_size(0),
	m_ndc(0), 
	m_debug(false),
	m_set_noise(false),   
	m_adc_ch(0), m_adc_xmin(0), m_adc_xmax(4000),
	m_prev_evtid(0),m_mv_ymin(-1000),m_mv_ymax(1000)
{
	gStyle->SetTitleFont(42,"Y");
	gStyle->SetTitleFont(42,"X");
	gStyle->SetLabelFont(42,"X");
	gStyle->SetLabelFont(42,"Y");
	gStyle->SetPalette(55);

	//
	m_dcevt_coll.clear();

	// Registration: InPort/OutPort/Service

	// Set InPort buffers
	registerInPort ("niutpcmonitor_in",  m_InPort);

	init_command_port();
	init_state_table();
	set_comp_name("NIUTPCMONITOR");

	// Initialize histgrams    
	for(int bo=0;bo<N_BOARD;bo++){
		for(int i=0; i<N_STRIP; i++) {
			// HG histogram
			m_h_adc_HG[bo][i] = new TH1D(Form("h_adc_bo%d_HGch%d",bo,i), Form("adc HGch%d",i), N_SAMPLE, 0, N_SAMPLE);
			m_h_adc_HG[bo][i]->SetStats(0);
			assert(m_h_adc_HG[bo][i]);
			
			m_h_adc_ped_HG[bo][i] = new TH1D(Form("h_adc_ped_bo%d_HGch%d",bo,i),Form("adc pedestal HGch%d",i), 4000, -1000,1000);
			assert(m_h_adc_ped_HG[bo][i]);


			// LG histogram
			m_h_adc_LG[bo][i] = new TH1D(Form("h_adc_bo%d_LGch%d",bo,i), Form("adc LGch%d",i), N_SAMPLE, 0, N_SAMPLE);
			m_h_adc_LG[bo][i]->SetStats(0);
			assert(m_h_adc_LG[bo][i]);

			m_h_adc_ped_LG[bo][i] = new TH1D(Form("h_adc_ped_bo%d_LGch%d",bo,i),Form("adc pedestal LGch%d",i), 4000, -1000,1000);
			assert(m_h_adc_ped_LG[bo][i]);
		}

		m_h_strips_HG[bo] = new TH2D(Form("m_h_strips_bo%d_HG",bo),"",N_SAMPLE,0,N_SAMPLE,N_STRIP,0,N_STRIP);
		m_h_strips_HG[bo]->SetStats(0);
		assert(m_h_strips_HG[bo]);

		m_h_pedestal_HG[bo] = new TH1D(Form("m_h_pedestal_bo%d_HG",bo),Form("pedestal BOARD:%d HG",bo),N_STRIP,0,N_STRIP);
		m_h_pedestal_HG[bo]->SetStats(0);
		m_h_pedestal_HG[bo]->GetYaxis()->SetTitle("[mV]");
		m_h_pedestal_HG[bo]->GetXaxis()->SetTitle("[CHANNEL]");      
		m_h_pedestal_HG[bo]->SetFillColor(kMagenta-2);
		assert(m_h_pedestal_HG[bo]);

		m_h_noise_HG[bo] = new TH1D(Form("m_h_noise_bo%d_HG",bo),Form("noise B0ARD:%d HG",bo),N_STRIP,0,N_STRIP);
		m_h_noise_HG[bo]->GetYaxis()->SetTitle("[mV]");
		m_h_noise_HG[bo]->GetXaxis()->SetTitle("[CHANNEL]");
		m_h_noise_HG[bo]->SetFillColor(kAzure-2);      
		m_h_noise_HG[bo]->SetStats(0);
		assert(m_h_noise_HG[bo]);

		m_h_strips_LG[bo] =
			new TH2D(Form("m_h_strips_bo%d_LG",bo),"",N_SAMPLE,0,N_SAMPLE,N_STRIP,0,N_STRIP);
		m_h_strips_LG[bo]->SetStats(0);
		assert(m_h_strips_LG[bo]);

		m_h_pedestal_LG[bo] =
			new TH1D(Form("m_h_pedestal_bo%d_LG",bo),Form("pedestal BOARD:%d LG",bo),N_STRIP,0,N_STRIP);
		m_h_pedestal_LG[bo]->SetStats(0);
		m_h_pedestal_LG[bo]->GetYaxis()->SetTitle("[mV]");
		m_h_pedestal_LG[bo]->GetXaxis()->SetTitle("[CHANNEL]");      
		m_h_pedestal_LG[bo]->SetFillColor(kMagenta-2);
		assert(m_h_pedestal_LG[bo]);

		m_h_noise_LG[bo] =
			new TH1D(Form("m_h_noise_bo%d_LG",bo),Form("noise B0ARD:%d LG",bo),N_STRIP,0,N_STRIP);
		m_h_noise_LG[bo]->GetYaxis()->SetTitle("[mV]");
		m_h_noise_LG[bo]->GetXaxis()->SetTitle("[CHANNEL]");
		m_h_noise_LG[bo]->SetFillColor(kAzure-2);      
		m_h_noise_LG[bo]->SetStats(0);
		assert(m_h_noise_LG[bo]);
	}

}

NIuTPCMonitor::~NIuTPCMonitor()
{
}


RTC::ReturnCode_t 
NIuTPCMonitor::onInitialize()
{
	if (m_debug) {
		std::cerr << "NIuTPCMonitor::onInitialize()" << std::endl;
	}

	return RTC::RTC_OK;
}

RTC::ReturnCode_t 
NIuTPCMonitor::onExecute(RTC::UniqueId ec_id)
{
	daq_do();

	return RTC::RTC_OK;
}

/*-----------------------daq_dummy------------------------------*/

int 
NIuTPCMonitor::daq_dummy()
{
	//
	int sleep_on = 0;
	for(int i=0;i<N_BOARD;i++){
		if(m_c[i]){
			// daq_dummy() will be invoked again after 10 msec.
			// This sleep reduces X servers' load.
			sleep_on = 1;
		}
	}

	if(sleep_on == 1){
		sleep(1);
		sleep_on = 0;
	}

	return 0;
}

/*-------------------------------------------------------------*/

/*-------------------daq_configure-----------------------------*/

int 
NIuTPCMonitor::daq_configure()
{
	std::cerr << "*** NIuTPCMonitor::configure" << std::endl;

	::NVList* paramList;
	paramList = m_daq_service0.getCompParams();
	parse_params(paramList);

	return 0;
}

/*-------------------------------------------------------------*/

/*********************parse_params*****************************/

int 
NIuTPCMonitor::parse_params(::NVList* list)
{

	std::cerr << "param list length:" << (*list).length() << std::endl;

	int len = (*list).length();
	for (int i = 0; i < len; i+=2) {
		std::string sname  = (std::string)(*list)[i].value;
		std::string svalue = (std::string)(*list)[i+1].value;

		std::cerr << "sname: " << sname << "  ";
		std::cerr << "value: " << svalue << std::endl;

		if (sname == "monitorUpdateRate") {
			if (m_debug) {
				std::cerr << "monitor update rate: " << svalue << std::endl;
			}
			char *offset;
			m_monitor_update_rate = (int)strtol(svalue.c_str(), &offset, 10);
		}
		// If you have more param in config.xml, write here
		else if (sname == "adcXmin")  {
			if (m_debug) { 
				std::cerr << "adc x min = " << svalue << std::endl;
			}
			char* offset;
			m_adc_xmin = (int)strtol(svalue.c_str(), &offset, 0);
		}
		else if (sname == "adcXmax")  {
			if (m_debug) { 
				std::cerr << "adc x max = " << svalue << std::endl;
			}
			char* offset;
			m_adc_xmax = (int)strtol(svalue.c_str(), &offset, 0);
		}
		else if (sname == "adcCh")  {
			if (m_debug) { 
				std::cerr << "adc ch = " << svalue << std::endl;
			}
			char* offset;
			m_adc_ch = (int)strtol(svalue.c_str(), &offset, 0);
		}
		else if (sname == "ndc")  {
			if (m_debug) { 
				std::cerr << "ndc = " << svalue << std::endl;
			}
			char* offset;
			m_ndc = (int)strtol(svalue.c_str(), &offset, 0);
		}
		else if (sname == "mv_ymin")  {
			if (m_debug) { 
				std::cerr << "mv_ymin = " << svalue << std::endl;
			}
			char* offset;
			m_mv_ymin = atof(svalue.c_str());
		}

		else if (sname == "mv_ymax")  {
			if (m_debug) { 
				std::cerr << "mv_ymax = " << svalue << std::endl;
			}
			char* offset;
			m_mv_ymax = atof(svalue.c_str());
		}
		else if (sname == "pedestal_window")  {
			if (m_debug) { 
				std::cerr << "pedestal_window = " << svalue << std::endl;
			}
			char* offset;
			m_pedestal_window = atoi(svalue.c_str());
		}
		else if (sname == "setNoise") {
			if (m_debug) {
				std::cerr << "setNoise: " << svalue << std::endl;
			}
			if (svalue == "yes" || svalue == "Yes" || svalue == "YES") {
				m_set_noise = true;
			}
			else {
				m_set_noise = false;
			}
		}
	}

	if (m_monitor_update_rate == 0) {
		m_monitor_update_rate = 1000;
	}

	std::cout << "Monitor parameters : " << std::endl
              << "  ndc          = " << m_ndc << std::endl
              << "  update_rate  = " << m_monitor_update_rate << std::endl;

	return 0;
}

/***************************************************************/

/*--------------------daq_uncofigure---------------------------*/

int 
NIuTPCMonitor::daq_unconfigure()
{

	std::cerr << "*** NIuTPCMonitor::unconfigure" << std::endl;
	//////////Delete Canvas//////////
	for(int i=0;i<N_BOARD;i++){
		if (m_c[i]) {
			delete m_c[i];
			m_c[i] = 0;
		}
	}
	if (m_c_noise){
		delete m_c_noise;
		m_c_noise =0;
	}
	return 0;
}

/*------------------------------------------------------------*/


/*------------------daq_start---------------------------------*/
int 
NIuTPCMonitor::daq_start()
{
	std::cerr << "*** NIuTPCMonitor::start" << std::endl;

	m_in_status  = BUF_SUCCESS;

	//////////Delete Canvas//////////
	for(int i=0;i<N_BOARD;i++){
		if (m_c[i]) {
			delete m_c[i];
			m_c[i] = 0;
		}
		m_c[i] = new TCanvas(Form("c_%d",i),"",1200,600);
		m_c[i]->Divide(2,3); // 4x2 plots		      
	}
	if( m_c_noise ){
		delete m_c_noise;
		m_c_noise=0;
	}
	if( m_set_noise ){
		m_c_noise = new TCanvas("c_noise","noise monitor",1200,600);
		m_c_noise->Divide(4,4);
	}

	return 0;
}

/*---------------------------------------------------------*/

/*--------------------daq_stop-----------------------------*/

int 
NIuTPCMonitor::daq_stop()
{
	std::cerr << "*** NIuTPCMonitor::stop" << std::endl;

	for(int i=0;i<N_BOARD;i++){
		if(m_c[i])
			m_c[i] -> Update();
	}
	if(m_c_noise && m_set_noise){
		m_c_noise->Update();
	}

	///////////////////
	reset_InPort();

	return 0;
}

/*---------------------------------------------------------*/

/*---------------------daq_pause---------------------------*/

int 
NIuTPCMonitor::daq_pause()
{
	std::cerr << "*** NIuTPCMonitor::pause" << std::endl;
	return 0;
}

/*---------------------------------------------------------*/

/*--------------------daq_resume---------------------------*/

int 
NIuTPCMonitor::daq_resume()
{
	std::cerr << "*** NIuTPCMonitor::resume" << std::endl;
	return 0;
}

/*---------------------------------------------------------*/

/**********************reset_InPort*************************/

int 
NIuTPCMonitor::reset_InPort()
{
	int ret = true;
	while(ret == true) {
		ret = m_InPort.read();
	}

	return 0;
}

/**********************************************************/

/*------------------decode_data---------------------------*/

int 
NIuTPCMonitor::decode_data(const unsigned char* mydata)
{
	return 0;
}

/*-----------------------------------------------------------*/

/******************************************************************/
/*---------------------user_main_logic----------------------------*/
/******************************************************************/

int 
NIuTPCMonitor::fill_data(const unsigned char* mydata, const int size)
{
	//
	assert(size%4==0);

	unsigned int* p = (unsigned int*)mydata;
	unsigned int* int_p = (unsigned int *)&mydata[8];
	int length = ntohl(*int_p);

	int data_set = length / 2 / N_CHANNEL;
	//std::cout << "data_set = " << data_set << std::endl;
	if( data_set > N_SAMPLE ) {
		std::cerr << "Too large # of sample" << std::endl;
		return 0;
	}

	//don't used
	dc_header header;
	header.m_magic    = p[0];
	header.m_id       = p[1];
	header.m_length   = p[2];
	header.m_trig_cnt = p[3];

	int mod_num=ntohl(p[1]);
	mod_num = mod_num & 0xffffff;

	std::cerr << "ID: "<< mod_num << "\t" << header.trig_cnt() << std::endl;    

	unsigned char* ptr_body = (unsigned char*)&(p[4]);
	m_dcevt_coll.push_back(dc_event(header, ptr_body, header.length()));

	int index=mod_num;

	// reset histo
	for(int i=0; i<N_STRIP; i++) {
		m_h_adc_HG[index][i]->Reset();
		m_h_adc_LG[index][i]->Reset();
		m_h_adc_ped_HG[index][i]->Reset();
		m_h_adc_ped_LG[index][i]->Reset();
	}
	m_h_strips_HG[index]->Reset();
	m_h_pedestal_HG[index]->Reset();
	m_h_noise_HG[index]->Reset();        
	m_h_strips_LG[index]->Reset();
	m_h_pedestal_LG[index]->Reset();
	m_h_noise_LG[index]->Reset();        


	std::vector< std::vector<double> > adc_data(N_CHANNEL,std::vector<double>(N_SAMPLE,0));

	////////////////////////////////////////
	// raw data format --> adc array
	for(int i = 0; i < data_set; i++){
		for(int ch = 0; ch < N_CHANNEL; ch++){
			int pos = 2*N_CHANNEL*i + 2*ch;
			unsigned short* short_p = (unsigned short *)&mydata[pos+16];
			adc_data[ch][i] = ntohs(*short_p);
		}
	}

	// pedestal calc
	double pedestals[N_CHANNEL];    
	for(int i=0;i<N_CHANNEL; i++){
		std::vector<double>::iterator itr=adc_data[i].begin();
		itr=itr+m_pedestal_window;
		pedestals[i]=(std::accumulate(adc_data[i].begin(),itr,0.0))/m_pedestal_window;
	}

	// adc array --> histo
	// 20171025 to mV (12bits , 2Vpp)
	int strip = 0;
	for(int ch=0; ch<N_CHANNEL; ch++) {
		for(int i=0; i<data_set;i++) {
			if(ch%2==0){
				m_h_adc_HG[index][strip]->SetBinContent( i+1, adc_data[ch][i]/pow(2,12)*2000-1000 );
				m_h_adc_ped_HG[index][strip]->Fill( adc_data[ch][i]/pow(2,12)*2000-1000 );
				double pulse_height=fabs((adc_data[ch][i]/pow(2,12)*2000-1000)-(pedestals[ch]/pow(2,12)*2000-1000));
				m_h_strips_HG[index]->Fill(i,strip,pulse_height );
			}
			if(ch%2==1){
				m_h_adc_LG[index][strip]->SetBinContent( i+1, adc_data[ch][i]/pow(2,12)*2000-1000 );
				m_h_adc_ped_LG[index][strip]->Fill( adc_data[ch][i]/pow(2,12)*2000-1000 );
				double pulse_height=fabs((adc_data[ch][i]/pow(2,12)*2000-1000)-(pedestals[ch]/pow(2,12)*2000-1000));
				m_h_strips_LG[index]->Fill(i,strip,pulse_height );
			}
		}
		if(ch%2==1)strip++;
	}

	//HG
	m_c[index]->cd(1);
	m_h_adc_HG[index][m_adc_ch]->GetXaxis()->SetRangeUser(m_adc_xmin, m_adc_xmax);
	m_h_adc_HG[index][m_adc_ch]->SetMinimum(m_mv_ymin);
	m_h_adc_HG[index][m_adc_ch]->SetMaximum(m_mv_ymax);
	m_h_adc_HG[index][m_adc_ch]->SetTitle(Form("Board: %d HG ch%d  Trigger:%d",mod_num,m_adc_ch,header.trig_cnt()));
	m_h_adc_HG[index][m_adc_ch]->GetXaxis()->SetTitle("sampling");
	m_h_adc_HG[index][m_adc_ch]->GetYaxis()->SetTitle("[mV]");     
	m_h_adc_HG[index][m_adc_ch]->DrawCopy();
	//
	m_c[index]->cd(3);
	m_h_strips_HG[index]->GetZaxis()->SetRangeUser(0,800);
	m_h_strips_HG[index]->SetContour(100);
	m_h_strips_HG[index]->Draw("colz");
	//
	m_c[index]->cd(5);
	m_h_adc_HG[index][m_adc_ch]->DrawCopy();     
	for(int i=0;i<N_STRIP;i++){
		if( i != m_adc_ch )
			m_h_adc_HG[index][i]->Draw("same");
	}
	// LG
	m_c[index]->cd(2);
	m_h_adc_LG[index][m_adc_ch]->GetXaxis()->SetRangeUser(m_adc_xmin, m_adc_xmax);
	m_h_adc_LG[index][m_adc_ch]->SetMinimum(m_mv_ymin);
	m_h_adc_LG[index][m_adc_ch]->SetMaximum(m_mv_ymax);
	m_h_adc_LG[index][m_adc_ch]->SetTitle(Form("Board: %d LG ch%d  Trigger:%d",mod_num,m_adc_ch,header.trig_cnt()));
	m_h_adc_LG[index][m_adc_ch]->GetXaxis()->SetTitle("sampling");
	m_h_adc_LG[index][m_adc_ch]->GetYaxis()->SetTitle("[mV]");     
	m_h_adc_LG[index][m_adc_ch]->DrawCopy();
	//
	m_c[index]->cd(4);
	m_h_strips_LG[index]->GetZaxis()->SetRangeUser(0,800);
	m_h_strips_LG[index]->SetContour(100);
	m_h_strips_LG[index]->Draw("colz");
	//
	m_c[index]->cd(6);
	m_h_adc_LG[index][m_adc_ch]->DrawCopy();     
	for(int i=0;i<N_STRIP;i++){
		if( i != m_adc_ch )
			m_h_adc_LG[index][i]->Draw("same");
	}
	//
	//m_c[index]->cd(4);
	m_c[index]->Update();    //Canvas update

	//========== noise =========//
	if( m_set_noise ){
		for(int i=0;i<N_STRIP;i++){
			// HG
			double rms_HG=m_h_adc_ped_HG[index][i]->GetRMS();
			double mean_HG=m_h_adc_ped_HG[index][i]->GetMean();
			m_h_pedestal_HG[index]->SetBinContent(i+1,mean_HG);
			m_h_noise_HG[index]->SetBinContent(i+1,rms_HG);	
			// LG
			double rms_LG=m_h_adc_ped_LG[index][i]->GetRMS();
			double mean_LG=m_h_adc_ped_LG[index][i]->GetMean();
			m_h_pedestal_LG[index]->SetBinContent(i+1,mean_LG);
			m_h_noise_LG[index]->SetBinContent(i+1,rms_LG);	
		}
		// HG pedestal
		m_c_noise->cd(index*4+1);
		m_h_pedestal_HG[index]->Draw();
		// HG noise
		m_c_noise->cd(index*4+2);
		m_h_noise_HG[index]->Draw();
		// LG pedestal
		m_c_noise->cd(index*4+3);
		m_h_pedestal_LG[index]->Draw();
		// LG noise
		m_c_noise->cd(index*4+4);
		m_h_noise_LG[index]->Draw();

		m_c_noise->Update();      
	}
	return 0;
}

/*------------------------------------------------------------*/

/************************read_InPort***************************/

unsigned int NIuTPCMonitor::read_InPort()
{
	/////////////// read data from InPort Buffer ///////////////
	unsigned int recv_byte_size = 0;
	bool ret = m_InPort.read();

	//////////////////// check read status /////////////////////
	if (ret == false) { // false: TIMEOUT or FATAL
		m_in_status = check_inPort_status(m_InPort);
		if (m_in_status == BUF_TIMEOUT) { // Buffer empty.
			if (check_trans_lock()) {     // Check if stop command has come.
				set_trans_unlock();       // Transit to CONFIGURE state.
            }
        }
        else if (m_in_status == BUF_FATAL) { // Fatal error
            fatal_error_report(INPORT_ERROR);
        }
    }
    else {
        recv_byte_size = m_in_data.data.length();
    }

    if (m_debug) {
        std::cerr << "m_in_data.data.length():" << recv_byte_size
                  << std::endl;
    }

    return recv_byte_size;
}

/******************************************************************/

/*-------------------daq_run--------------------------------------*/

int NIuTPCMonitor::daq_run()
{
    if (m_debug) {
        std::cerr << "*** NIuTPCMonitor::run" << std::endl;
    }


    // read from in-port
    unsigned int recv_byte_size = read_InPort();
    if (recv_byte_size == 0) { // Timeout
        return 0;
    }

    check_header_footer(m_in_data, recv_byte_size); // check header and footer
    m_event_byte_size = get_event_size(recv_byte_size);
    if (m_event_byte_size > DATA_BUF_SIZE) {
        fatal_error_report(USER_DEFINED_ERROR1, "DATA BUF TOO SHORT");
    }

    /////////////  Write component main logic here. /////////////
    memcpy(&m_recv_data[0], &m_in_data.data[HEADER_BYTE_SIZE], m_event_byte_size);

    fill_data(&m_recv_data[0], m_event_byte_size);

    unsigned long sequence_num = get_sequence_num();
    //
    if( (int)m_dcevt_coll.size() == m_ndc ) {
        int evt_num = sequence_num%10;
        if ((evt_num % m_monitor_update_rate) == 0) {
            do_analysis();
        }
        m_dcevt_coll.clear();
    }

    /////////////////////////////////////////////////////////////
    inc_total_data_size(m_event_byte_size);  // increase total data byte size
    inc_sequence_num();                      // increase sequence num.


    return 0;
}

void 
NIuTPCMonitor::do_analysis()
{
    //
    std::cout << "do_analysis() seq# = " << get_sequence_num() << std::endl;
    assert(m_dcevt_coll.size()<=10);
    int ndata[10] = {0};
    int trignum[10] = {0};
    for(int i=0; i<(int)m_dcevt_coll.size(); i++) {
        int id = m_dcevt_coll[i].header().id()&0xffff;
        int trgcnt = m_dcevt_coll[i].header().trig_cnt();
        m_dcevt_coll[i].header().show();

        std::cout << "id = " << id << " trgcnt = " << trgcnt << std::endl;
	
    }
    std::cout << "================= end of event" << std::endl;
}

extern "C"
{
    void NIuTPCMonitorInit(RTC::Manager* manager)
    {
        RTC::Properties profile(niutpcmonitor_spec);
        manager->registerFactory(profile,
                                 RTC::Create<NIuTPCMonitor>,
                                 RTC::Delete<NIuTPCMonitor>);
    }
};
