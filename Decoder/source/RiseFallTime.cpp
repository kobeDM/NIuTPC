//------------------------------------------------------------------
// Read tree for NIuTPc output data
// Version 0.1
// Update: 29. July 2016
// Author: T.Ikeda
//------------------------------------------------------------------

// STL
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
//
#include <time.h>
// ROOT
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TClassTable.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
//USER
#include "Event.h"

#define N_CHANNEL 64
#define N_STRIP 32
#define ADC_RANGE 4096
#define N_SAMPLING 4000 //Hz
#define SAMPLING_HELZ 2.5e6 //Hz

#define DEBUG 0
using namespace std;

//------------- main -------------------------------------------------
int main(int argc,char *argv[]){
	clock_t t1,t2;

	if(argc <3){
		std::cerr << "Usage:" << std::endl;
		std::cerr << "/.main [filename.root] [sample window]" << std::endl;
		return 1;
	}

	string filename=argv[1];
	std::string::size_type index = filename.find(".root");
	if( index == std::string::npos ) { 
		std::cout << "Failure!!!" << std::endl;
		return 1;
	}

	//json database
	string conffilename=argv[2];
	std::string::size_type index_conf = conffilename.find(".json");
	if( index_conf == std::string::npos ) { 
		std::cout << "Failure!!!" << std::endl;
		return 1;
	}
	boost::property_tree::ptree pt;
	read_json(conffilename,pt);
	int offset_samp = pt.get<int>("config.offset_sampling");


	t1=clock();
	std::cerr << "======================================" << std::endl;
	std::cerr << "Read ROOT file" << std::endl;
	std::cerr << "Visualizer for 0.1c NIuTPC" << std::endl;
	std::cerr << "Version 0.1" << std::endl;
	std::cerr << "======================================" << std::endl;
	std::cerr << "file name: " << filename << endl;
	std::cerr << "config file name: " << conffilename << endl;
	std::cerr << "======================================" << std::endl;
	std::cerr << "offset sampling: " << offset_samp << endl;

	TApplication app("app",&argc,argv);  

	if(!TClassTable::GetDict("Event.so")){
		gSystem->Load("/home/msgc/LTARS2018/analyzer/Decoder0.2_cmake/bin/libEvent.so");
	}

	TFile *f=new TFile(filename.c_str());
	if(!f){
		std::cerr << "ERROR: Cant find file" << std::endl;
		return 1;
	}

	TTree *tree=(TTree*)f->Get("Tree");

	Event *event = new Event();
	int nevent=tree->GetEntries();
	std::cerr << "event number is " << nevent << std::endl;

	TBranch *branch = tree->GetBranch("Event");
	branch->SetAddress(&event);

	TH1F *h_hgtime[N_STRIP];
	TH1F *h_lgtime[N_STRIP];

	TH1F *h_hgrise[N_STRIP];
	TH1F *h_lgrise[N_STRIP];
	
	TH1F *h_hgfall[N_STRIP];
	TH1F *h_lgfall[N_STRIP];

	for(int st=0;st<N_STRIP;st++){
		h_hgtime[st] = new TH1F(Form("h_hgtimeCH%d",st),    Form("HighGainPeakTimeCH%d",st),N_SAMPLING,0,N_SAMPLING);
		h_lgtime[st] = new TH1F(Form("h_lgtimeCH%d",st),    Form("LowGainPeakTimeCH%d", st),N_SAMPLING,0,N_SAMPLING);
		h_hgrise[st] = new TH1F(Form("h_hgriseCH%d",st),    Form("HighGainRiseTimeCH%d",st),N_SAMPLING,0,N_SAMPLING);
		h_lgrise[st] = new TH1F(Form("h_lgriseCH%d",st),    Form("LowGainRiseTimeCH%d", st),N_SAMPLING,0,N_SAMPLING);
		h_hgfall[st] = new TH1F(Form("h_hgfallCH%d",st),    Form("HighGainFallTimeCH%d",st),N_SAMPLING,0,N_SAMPLING);
		h_lgfall[st] = new TH1F(Form("h_lgfallCH%d",st),    Form("LowGainFallTimeCH%d", st),N_SAMPLING,0,N_SAMPLING);
	}

	//---------------EventLoop--------------------
	for(int ev=0;ev<nevent;ev++){
	//for(int ev=0;ev<1;ev++){

		tree->GetEntry(ev);

		int event_id               = event->GetEventID();
		int sampling_num           = event->GetSamplingNum();
		vector<vector<double>> adc = event->GetADC();

		//---------------Waveform--------------------
		//Even CH -> High Gain
		//Odd CH  -> Low  Gain
		
		int strip_num = 0;
		for(int j=0;j<N_CHANNEL;j++){
			double sum  = 0;
			double offset = 0;
			double height = 0;
			double time = 0;
			double rise_time_on = 0;
			double rise_time_off = 0;
			double fall_time = 0;
			bool rise_on = false;
			bool rise_off = false;
			bool fall_on = false;

			for(int k=0;k<offset_samp;k++){
				//HG
				if(j%2==0){
					sum += adc.at(j).at(k);
				}
				//LG
				if(j%2==1){
					sum += adc.at(j).at(k);
				}
			}
	
			offset = sum / double(offset_samp);
			std::vector<double>::iterator itr_max = std::max_element(adc.at(j).begin(),adc.at(j).end());
			size_t index = std::distance(adc.at(j).begin(), itr_max);
			height = adc.at(j)[index]-offset;
			time = double(index);
			for(int k=0;k<sampling_num;k++){
				if(k > sampling_num-index && adc.at(j).at(sampling_num-k-1)-offset <= height*0.1){
					if(!rise_on){
						rise_on=true;
						rise_time_on = k;
					}
				}
				if(k > sampling_num-index && adc.at(j).at(sampling_num-k-1)-offset <= height*0.9){
					if(!rise_off){
						rise_off=true;
						rise_time_off = k;
					}
				}
				if(k>index && adc.at(j).at(k)-offset <= height/std::exp(1)){
					if(!fall_on){
						fall_on=true;
						fall_time = k;
					}
				}
			}

			//HG
			if(j%2==0){
				h_hgtime[strip_num]->Fill(time);
				h_hgrise[strip_num]->Fill(rise_time_on-rise_time_off);
				h_hgfall[strip_num]->Fill(fall_time-time);
			}
			//LG
			if(j%2==1){
				h_lgtime[strip_num]->Fill(time);
				h_lgrise[strip_num]->Fill(rise_time_on-rise_time_off);
				h_lgfall[strip_num]->Fill(fall_time-time);
				strip_num++;
			}
		}
	}

	TGraphErrors* ge_hg_rise = new TGraphErrors();
	ge_hg_rise->SetTitle("HighGainCHRiseTime");
	ge_hg_rise->GetXaxis()->SetTitle("ch");
	ge_hg_rise->GetYaxis()->SetTitle("Clock(4000sampling/2.5MHz)");
	TGraphErrors* ge_lg_rise = new TGraphErrors();
	ge_lg_rise->SetTitle("LowGainCHRiseTime");
	ge_lg_rise->GetXaxis()->SetTitle("ch");
	ge_lg_rise->GetYaxis()->SetTitle("Clock(4000sampling/2.5MHz)");
	
	TGraphErrors* ge_hg_fall = new TGraphErrors();
	ge_hg_fall->SetTitle("HighGainCHFallTime");
	ge_hg_fall->GetXaxis()->SetTitle("ch");
	ge_hg_fall->GetYaxis()->SetTitle("Clock(4000sampling/2.5MHz)");
	TGraphErrors* ge_lg_fall = new TGraphErrors();
	ge_lg_fall->SetTitle("LowGainCHFallTime");
	ge_lg_fall->GetXaxis()->SetTitle("ch");
	ge_lg_fall->GetYaxis()->SetTitle("Clock(4000sampling/2.5MHz)");

	TGraphErrors* ge_hg_time = new TGraphErrors();
	ge_hg_time->SetTitle("HighGainCHtime");
	ge_hg_time->GetXaxis()->SetTitle("ch");
	ge_hg_time->GetYaxis()->SetTitle("Clock(4000sampling/2.5MHz)");
	TGraphErrors* ge_lg_time = new TGraphErrors();
	ge_lg_time->SetTitle("LowGainCHtime");
	ge_lg_time->GetXaxis()->SetTitle("ch");
	ge_lg_time->GetYaxis()->SetTitle("Clock(4000sampling/2.5MHz)");
	
	//---------------Draw--------------------
	//TCanvas *c_strip=new TCanvas("c_strip","",800,800);
	//h_strip->Draw("colz");
	
	//PeakTime
	TCanvas *c_pt=new TCanvas("c_pt","c_pt",1000,1000);
	c_pt->Divide(8,8);
	for(int st=0;st<N_STRIP;st++){

		h_hgtime[st]->GetXaxis()->SetTitle("Clock (4000sampling/2.5MHz)");
		h_hgtime[st]->GetYaxis()->SetTitle("count");
		h_lgtime[st]->GetXaxis()->SetTitle("Clock (4000sampling/2.5MHz)");
		h_lgtime[st]->GetYaxis()->SetTitle("count");
		
		c_pt->cd(st+1);    h_hgtime[st]->Draw("hist");
		c_pt->cd(st+32+1); h_lgtime[st]->Draw("hist");
		
		ge_hg_time->SetPoint(     st, st, h_hgtime[st]->GetMean());
		ge_hg_time->SetPointError(st, 0,  h_hgtime[st]->GetMeanError());
		ge_lg_time->SetPoint(     st, st, h_lgtime[st]->GetMean());
		ge_lg_time->SetPointError(st, 0,  h_lgtime[st]->GetMeanError());
	}

	TCanvas *c_time = new TCanvas("c_time","c_time",1000,500);
	c_time->Divide(2,1);
	c_time->cd(1);
	ge_hg_time->Draw("ap");
	c_time->cd(2);
	ge_lg_time->Draw("ap");

	//RiseTime
	TCanvas *c_rise=new TCanvas("c_rise","c_rise",1000,1000);
	c_rise->Divide(8,8);
	for(int st=0;st<N_STRIP;st++){

		h_hgrise[st]->GetXaxis()->SetTitle("Clock (4000sampling/2.5MHz)");
		h_hgrise[st]->GetYaxis()->SetTitle("count");
		h_lgrise[st]->GetXaxis()->SetTitle("Clock (4000sampling/2.5MHz)");
		h_lgrise[st]->GetYaxis()->SetTitle("count");
		
		c_rise->cd(st+1);    h_hgrise[st]->Draw("hist");
		c_rise->cd(st+32+1); h_lgrise[st]->Draw("hist");
		
		ge_hg_rise->SetPoint(     st, st, h_hgrise[st]->GetMean());
		ge_hg_rise->SetPointError(st, 0,  h_hgrise[st]->GetMeanError());
		ge_lg_rise->SetPoint(     st, st, h_lgrise[st]->GetMean());
		ge_lg_rise->SetPointError(st, 0,  h_lgrise[st]->GetMeanError());
	}

	TCanvas *c_RISE = new TCanvas("c_RISE","c_RISE",1000,500);
	c_RISE->Divide(2,1);
	c_RISE->cd(1);
	ge_hg_rise->Draw("ap");
	c_RISE->cd(2);
	ge_lg_rise->Draw("ap");

	//FallTime
	TCanvas *c_fall=new TCanvas("c_fall","c_fall",1000,1000);
	c_fall->Divide(8,8);
	for(int st=0;st<N_STRIP;st++){

		h_hgfall[st]->GetXaxis()->SetTitle("Clock (4000sampling/2.5MHz)");
		h_hgfall[st]->GetYaxis()->SetTitle("count");
		h_lgfall[st]->GetXaxis()->SetTitle("Clock (4000sampling/2.5MHz)");
		h_lgfall[st]->GetYaxis()->SetTitle("count");
		
		c_fall->cd(st+1);    h_hgfall[st]->Draw("hist");
		c_fall->cd(st+32+1); h_lgfall[st]->Draw("hist");
		
		ge_hg_fall->SetPoint(     st, st, h_hgfall[st]->GetMean());
		ge_hg_fall->SetPointError(st, 0,  h_hgfall[st]->GetMeanError());
		ge_lg_fall->SetPoint(     st, st, h_lgfall[st]->GetMean());
		ge_lg_fall->SetPointError(st, 0,  h_lgfall[st]->GetMeanError());
	}

	TCanvas *c_FALL = new TCanvas("c_FALL","c_FALL",1000,500);
	c_FALL->Divide(2,1);
	c_FALL->cd(1);
	ge_hg_fall->Draw("ap");
	c_FALL->cd(2);
	ge_lg_fall->Draw("ap");

	//---------------Write--------------------
	std::string outfilename = filename.substr(0,index) + "_time.root";
	TFile* outfile = new TFile(outfilename.c_str(),"recreate");

	for(int st=0;st<N_STRIP;st++){
		h_hgtime[st]->Write();
		h_lgtime[st]->Write();
		h_hgrise[st]->Write();
		h_lgrise[st]->Write();
		h_hgfall[st]->Write();
		h_lgfall[st]->Write();
	}
	ge_hg_time->Write();
	ge_lg_time->Write();
	ge_hg_rise->Write();
	ge_lg_rise->Write();
	ge_hg_fall->Write();
	ge_lg_fall->Write();

	outfile->Close();


	//----------- End of output ------------//
	std::cerr << std::endl;
	std::cerr << "======================================" << std::endl;
	t2=clock();
	std::cerr << "execution time: " << (t2-t1)/CLOCKS_PER_SEC << std::endl;

	app.Run();

	return 1;
}

