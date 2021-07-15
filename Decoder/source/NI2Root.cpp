//============================================================
// Decorder for NIuTPC Logger
//
// Input data format is 
//   HEADER(16byte) + ADCDATA(4000*2*64) + TIMESTAMP(4byte)
//------------------------------------------------------------
// Update : 29. July 2016
// Author : T.Ikeda
//============================================================


#include <stdio.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//ROOT
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TApplication.h"
// USER
#include "Event.h"

// parameters
#define HEADER_SIZE 16 //byte
#define N_CHANNEL 64
#define N_BOARD_STRIP 32
#define N_SAMPLE 4000  // always 4000
#define SAMPLING_HELZ 2.5e6 //default

using namespace std;

struct event_info {
	int module_num = -1;
	int trigger;
	int timestamp;
	vector<vector<double>> hg_adc;
	vector<vector<double>> lg_adc;
};

//#########################
// main
//#########################
int main(int argc,char *argv[]){
	
	gStyle->SetPalette(kRainBow);

	std::cerr << "=========================================" << std::endl;
	std::cerr << " Create ROOT file for NIuTPC output data " << std::endl;
	std::cerr << " Versioin 0.1 " << std::endl;
	std::cerr << " Using : ./main [filename] " << std::endl;
	std::cerr << "=========================================" << std::endl;

	if(argc != 2) {
		std::cout << "<ERROR>    Specify input file name as first argument." 
                  << std::endl;
		return 1;
	}

	string dirfilename=argv[1];
	std::string::size_type index = dirfilename.find(".dat");
	if( index == std::string::npos ) { 
		std::cout << "Failure!!!" << std::endl;
		return 1;
	}

	// ------- get file size ------//
	FILE *fp;
	long file_size;
	char *buffer;
	struct stat stbuf;
	int fd;
	fd = open(dirfilename.c_str(), O_RDONLY);
	if (fd == -1) {
		cerr << "ERROR: Cannot get file size" << endl;
		return -1;
	}
	fp = fdopen(fd, "rb");
	if (fp == NULL) {
		cerr << "ERROR: Cannot get file size" << endl;
		return -1;
	}
	if (fstat(fd, &stbuf) == -1) {
		cerr << "ERROR: Cannot get file size" << endl;
		return -1;
	}
	file_size = stbuf.st_size;
	cerr << "Data size: " << file_size << "byte" << endl;


	//open file
	ifstream fin( dirfilename.c_str(), ios::in | ios::binary );
	if (!fin){
		cerr << "ERROR: Cannot open file" << endl;
		return 1;
	}
	
	//----------- Creat root file --------//
	std::string filename;
	std::string::size_type pos = dirfilename.find("/");
	while(pos != std::string::npos){
		dirfilename = dirfilename.substr(pos+1);
		pos = dirfilename.find("/");
	}
	pos = dirfilename.find(".dat");
	filename = dirfilename.substr(0,pos);
	std::string outfilename = filename.substr(0, index) + ".root";
	TFile *f=new TFile(outfilename.c_str(),"recreate");
	// ---------- Create Tree ------------//
	TTree* tree= new TTree("Tree","Event Tree");

	// ---------- Create Event Class -----//
	Event *event = new Event();
	tree->Branch("Event", "Event", &event);

	vector<event_info> anode_info_0(3000);
	vector<event_info> anode_info_1(3000);
	vector<event_info> cathode_info_0(3000);
	vector<event_info> cathode_info_1(3000);

	/*/monitor
      TApplication app("app",&argc,argv);  
      TCanvas* c = new TCanvas("c","c",1200,800);
      c->Divide(2,2);
      c->cd(1)->DrawFrame(0,-1000,4000,1000,"waveform HG;sampling;mV");
      c->cd(2)->DrawFrame(0,-1000,4000,1000,"waveform LG;sampling;mV");
	*/
	
	double current_size=0;
	int ev_num=0;

	int first_trigger_num;

	while(!fin.eof()){
    
		if(ev_num%1==0) std::cout << "\rFill Data into BUFFER ... : " 
                                  << ev_num << "/" << double(file_size/512020) << std::flush;
		// header read
		unsigned char header[HEADER_SIZE];
		fin.read( ( char * ) &header[0], HEADER_SIZE*sizeof( char ) );
        unsigned int* p = (unsigned int *)&header;
        int m_magic = p[0];
        // ID
        int module_num = ntohl(p[1]);
        module_num = module_num & 0xffffff;
        // DATA length
        int length = ntohl(p[2]); // unit?
        // trigger count
        int trigger = ntohl(p[3]);
        int data_set = length/2/N_CHANNEL;
        if(data_set > N_SAMPLE){
            cerr << "ERROR: Too large # of sample" << endl;
            return 1;
        }
        //data read
        char *data_s;
        data_s=new char[length];
        fin.read( ( char * ) &data_s[0], length*sizeof(char ) );

        // time stamp read
        int timestamp;
        fin.read( ( char * ) &timestamp, sizeof(int ) );
        timestamp = ntohl(timestamp);

		int board_strip=0;
		vector< vector<double> > tmp_hg_adc(N_BOARD_STRIP,vector<double>(N_SAMPLE,0));
		vector< vector<double> > tmp_lg_adc(N_BOARD_STRIP,vector<double>(N_SAMPLE,0));
		for(int ch = 0;ch < N_CHANNEL;ch++){
			for(int i=0;i<data_set;i++){
				int pos = 2*N_CHANNEL*i + 2*ch;
				unsigned short* short_p = (unsigned short *)&data_s[pos];
				if(ch%2==0){
					tmp_hg_adc[board_strip][i] = ntohs(*short_p)/pow(2,12)*2000-1000; // mV
				}
				if(ch%2==1){
					tmp_lg_adc[board_strip][i] = ntohs(*short_p)/pow(2,12)*2000-1000; // mV
				}
			}
			if(ch%2==1)board_strip++;
		}
  
        /*
          cout << " Trigger : " << trigger << " Module     : " << module_num
          << " Length  : " << length  << " Time Stamp : " << timestamp
          << " 1st HG mV : " << tmp_hg_adc[0][0]
          << " 1st LG mV : " << tmp_lg_adc[0][0]
          << endl;
        */
		/*/ monitor
          TGraph* g_wf_hg[32];
          TGraph* g_wf_lg[32];
          for(int i=0;i<32;i++){
          g_wf_hg[i] = new TGraph();
          g_wf_lg[i] = new TGraph();
          }
          TH2D*   h_strip_hg = new TH2D("h_strip_hg","h_strip_hg",4000,0,4000,32,0,32);
          TH2D*   h_strip_lg = new TH2D("h_strip_lg","h_strip_lg",4000,0,4000,32,0,32);
          board_strip=0;
          for(int ch = 0;ch < N_CHANNEL;ch++){
          for(int i=0;i<data_set;i++){
          if(ch%2==0){
          g_wf_hg[board_strip]->SetPoint(i,i,tmp_hg_adc[board_strip][i]);
          h_strip_hg->SetBinContent(i,board_strip,tmp_hg_adc[board_strip][i]);
          }
          if(ch%2==1){
          g_wf_lg[board_strip]->SetPoint(i,i,tmp_lg_adc[board_strip][i]);
          h_strip_lg->SetBinContent(i,board_strip,tmp_lg_adc[board_strip][i]);
          }
          }
          if(ch%2==1)board_strip++;
          }

          for(int i=0;i<32;i++){
          c->cd(1);g_wf_hg[i]->Draw("p same");
          c->cd(2);g_wf_lg[i]->Draw("p same");
          }
          c->cd(3);
          h_strip_hg->Draw("colz");
          c->cd(4);
          h_strip_lg->Draw("colz");

          if(trigger==1&&module_num==3)app.Run();
          for(int i=0;i<32;i++){
          delete g_wf_hg[i];
          delete g_wf_lg[i];
          }
          delete h_strip_hg;
          delete h_strip_lg;
		*/
		if(ev_num==0)first_trigger_num = trigger;

		// buffer fill
		event_info ev_tmp;
		ev_tmp.module_num = module_num;
		ev_tmp.timestamp = timestamp;
		ev_tmp.trigger = trigger;
		ev_tmp.hg_adc = tmp_hg_adc;
		ev_tmp.lg_adc = tmp_lg_adc;
		if(module_num==0){
			cathode_info_0[trigger-first_trigger_num+10] = ev_tmp;
		}
		else if(module_num==1){
			cathode_info_1[trigger-first_trigger_num+10] = ev_tmp;
		}
		else if(module_num==2){
			anode_info_0[trigger-first_trigger_num+10]   = ev_tmp;
		}
		else if(module_num==3){
			anode_info_1[trigger-first_trigger_num+10]   = ev_tmp;
		}
		current_size += double(HEADER_SIZE + length + 4);
		ev_num++;
	}
	
	double nevent = double(ev_num/4);
	ev_num=0;
	std::cout<<std::endl;

	for(int i=0;i<anode_info_0.size();i++){
		if(anode_info_0[i].module_num==-1 || 
           anode_info_1[i].module_num==-1 || 
           cathode_info_0[i].module_num==-1 || 
           cathode_info_1[i].module_num==-1 
           ){
			//std::cout << "skip trigger ev" <<std::endl;
			continue;
		}
		if(ev_num%1==0) std::cerr << "\rFill Data into Tree ... : "<< ev_num << "/" << nevent <<std::flush;
		ev_num++;
		//12.11 module2
		vector<vector<double>> hg_a_adc_0 = anode_info_0[i].hg_adc;
		vector<vector<double>> lg_a_adc_0 = anode_info_0[i].lg_adc;
		//13.11 module2
		vector<vector<double>> hg_a_adc_1 = anode_info_1[i].hg_adc;
		vector<vector<double>> lg_a_adc_1 = anode_info_1[i].lg_adc;
		//10.11 module2
		vector<vector<double>> hg_c_adc_0 = cathode_info_0[i].hg_adc;
		vector<vector<double>> lg_c_adc_0 = cathode_info_0[i].lg_adc;
		//11.11 module2
		vector<vector<double>> hg_c_adc_1 = cathode_info_1[i].hg_adc;
		vector<vector<double>> lg_c_adc_1 = cathode_info_1[i].lg_adc;
		
		//combine data
		/*
          hg_a_adc_0.reserve(hg_a_adc_0.size(),hg_a_adc_1.size());
          lg_a_adc_0.reserve(lg_a_adc_0.size(),lg_a_adc_1.size());
          hg_c_adc_0.reserve(hg_c_adc_0.size(),hg_c_adc_1.size());
          lg_c_adc_0.reserve(lg_c_adc_0.size(),lg_c_adc_1.size());
          std::copy(hg_a_adc_1.begin(), hg_a_adc_1.end(), std::back_inserter(hg_a_adc_0));
          std::copy(lg_a_adc_1.begin(), lg_a_adc_1.end(), std::back_inserter(lg_a_adc_0));
          std::copy(hg_c_adc_1.begin(), hg_c_adc_1.end(), std::back_inserter(hg_c_adc_0));
          std::copy(lg_c_adc_1.begin(), lg_c_adc_1.end(), std::back_inserter(lg_c_adc_0));
		*/
		hg_a_adc_0.insert(hg_a_adc_0.end(),hg_a_adc_1.begin(),hg_a_adc_1.end());
		lg_a_adc_0.insert(lg_a_adc_0.end(),lg_a_adc_1.begin(),lg_a_adc_1.end());
		hg_c_adc_0.insert(hg_c_adc_0.end(),hg_c_adc_1.begin(),hg_c_adc_1.end());
		lg_c_adc_0.insert(lg_c_adc_0.end(),lg_c_adc_1.begin(),lg_c_adc_1.end());

		event->SetHeaders(anode_info_0[i].trigger,anode_info_0[i].timestamp);
		event->SetAnodeHGADC(hg_a_adc_0);
		event->SetAnodeLGADC(lg_a_adc_0);
		event->SetCathodeHGADC(hg_c_adc_0);
		event->SetCathodeLGADC(lg_c_adc_0);

		tree->Fill();
	}
	
	tree->Write();
	f->Close(); delete f;
	cerr << endl;
	fin.close();
	return 0;

}
