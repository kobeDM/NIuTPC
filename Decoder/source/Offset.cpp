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


	TH1F *h_hgoff[N_STRIP];
	TH1F *h_lgoff[N_STRIP];

	for(int st=0;st<N_STRIP;st++){
		h_hgoff[st] = new TH1F(Form("h_hgoffCH%d",st),Form("HighGainOffsetCH%d",st),ADC_RANGE,0,ADC_RANGE);
		h_lgoff[st] = new TH1F(Form("h_lgoffCH%d",st),Form("LowGainOffsetCH%d",st) ,ADC_RANGE,0,ADC_RANGE);
	}

	//---------------EventLoop--------------------
	for(int ev=0;ev<nevent;ev++){

		tree->GetEntry(ev);

		int event_id               = event->GetEventID();
		int sampling_num           = event->GetSamplingNum();
		vector<vector<double>> adc = event->GetADC();

		//---------------Waveform--------------------
		//Even CH -> High Gain
		//Odd CH  -> Low  Gain
		
		int strip_num = 0;
		for(int j=0;j<N_CHANNEL;j++){
			//printf("debug%d\n",j);
			double sum  = 0;
			double mean = 0;
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
			mean = sum / double(offset_samp);
			if(j%2==0){
				h_hgoff[strip_num]->Fill(mean);
			}
			//LG
			if(j%2==1){
				h_lgoff[strip_num]->Fill(mean);
				strip_num++;
			}
		}
	}
	TGraphErrors* ge_hg_off = new TGraphErrors();
	ge_hg_off->SetTitle("HighGainCHOffset");
	ge_hg_off->GetXaxis()->SetTitle("ch");
	ge_hg_off->GetYaxis()->SetTitle("ADC");
	TGraphErrors* ge_lg_off = new TGraphErrors();
	ge_lg_off->SetTitle("LowGainCHOffset");
	ge_lg_off->GetXaxis()->SetTitle("ch");
	ge_lg_off->GetYaxis()->SetTitle("ADC");

	//---------------Draw--------------------
	//TCanvas *c_strip=new TCanvas("c_strip","",800,800);
	//h_strip->Draw("colz");
	TCanvas *c_off=new TCanvas("c_off","c_off",1000,1000);
	c_off->Divide(8,8);
	for(int st=0;st<N_STRIP;st++){

		h_hgoff[st]->GetXaxis()->SetTitle("ADC (-1~1V/12bit)");
		h_hgoff[st]->GetYaxis()->SetTitle("count");
		h_lgoff[st]->GetXaxis()->SetTitle("ADC (-1~1V/12bit)");
		h_lgoff[st]->GetYaxis()->SetTitle("count");
		
		c_off->cd(st+1);    h_hgoff[st]->Draw("hist");
		c_off->cd(st+32+1);	h_lgoff[st]->Draw("hist");
		
		ge_hg_off->SetPoint(     st, st, h_hgoff[st]->GetMean());
		ge_hg_off->SetPointError(st, 0,  h_hgoff[st]->GetMeanError());
		ge_lg_off->SetPoint(     st, st, h_lgoff[st]->GetMean());
		ge_lg_off->SetPointError(st, 0,  h_lgoff[st]->GetMeanError());
	}

	TCanvas *c_offset = new TCanvas("c_offset","c_offset",1000,500);
	c_offset->Divide(2,1);
	c_offset->cd(1);
	ge_hg_off->Draw("ap");
	c_offset->cd(2);
	ge_lg_off->Draw("ap");



	//---------------Write--------------------
	std::string outfilename = filename.substr(0,index) + "_off.root";
	TFile* outfile = new TFile(outfilename.c_str(),"recreate");

	for(int st=0;st<N_STRIP;st++){
		h_hgoff[st]->Write();
		h_lgoff[st]->Write();
	}
	ge_hg_off->Write();
	ge_lg_off->Write();

	outfile->Close();


	//----------- End of output ------------//
	std::cerr << std::endl;
	std::cerr << "======================================" << std::endl;
	t2=clock();
	std::cerr << "execution time: " << (t2-t1)/CLOCKS_PER_SEC << std::endl;

	app.Run();

	return 1;
}

