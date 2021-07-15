//------------------------------------------------------------------
// Read tree for NIuTPc output data
// Version 0.1
// Update: 27. August 2020
// Author: T.Shimada
//------------------------------------------------------------------

// STL
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
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
#include "TLegend.h"
#include "TSpectrum.h"
//USER
#include "Event.h"
#include "NIConfig.hh"

#define N_CHANNEL 64
#define N_STRIP 32
#define SAMPLING_NUM 4000
#define SAMPLING_HELZ 2.5e6 //Hz

#define DEBUG 0
using namespace std;

//------------- main -------------------------------------------------
int main(int argc,char *argv[]){
	clock_t t1,t2;
	

	Int_t MyPalette[64];
	Double_t Red   [] = {1.0, 0.0, 1.0};
	Double_t Green [] = {0.0, 0.0, 0.0};
	Double_t Blue  [] = {0.0, 1.0, 1.0};
	Double_t Length[] = {0.0, 0.5, 1.0};
	Int_t FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 64);
	for (int i=0;i<64;i++) MyPalette[i] = FI+i;

	if(argc <4){
		std::cerr << "Usage:" << std::endl;
		std::cerr << "/.main [filename.root] [config.json] [event number]" << std::endl;
		return 1;
	}
	
	//++++++++++++++++++++++++++++++++++++++++++
	//  read root file 
	//++++++++++++++++++++++++++++++++++++++++++
	std::string dirfilename=argv[1];
	std::string::size_type index = dirfilename.find(".root");
	if( index == std::string::npos ) { 
		std::cout << "Failure!!!" << std::endl;
		return 1;
	}
	//++++++++++++++++++++++++++++++++++++++++++
	//  read config file 
	//++++++++++++++++++++++++++++++++++++++++++
	string conffilename=argv[2];
	NIConfig* ni_conf = new NIConfig();
	ni_conf->ReadConfigJSON(conffilename);
	int    offset_sampling       = ni_conf->offset_sampling;
	double driftV_main           = ni_conf->driftV_main;
	double driftV_mino           = ni_conf->driftV_mino;
	double calc_abs_z_param      = driftV_main*driftV_mino/(driftV_mino-driftV_main);
	double tot_anode_threshold   = ni_conf->tot_anode_threshold;
	double lg_anode_threshold    = ni_conf->lg_anode_threshold;
	double tot_cathode_threshold = ni_conf->tot_cathode_threshold;
	double lg_cathode_threshold  = ni_conf->lg_cathode_threshold;
	double minority_threshold    = ni_conf->minority_threshold;
	double minority_ROI_start    = ni_conf->minority_ROI_start;
	double minority_ROI_end      = ni_conf->minority_ROI_end;

	//event number
	int ev_num = atoi(argv[3]);

	//++++++++++++++++++++++++++++++++++++++++++
	//  analysis start
	//++++++++++++++++++++++++++++++++++++++++++
	std::cerr << "======================================" << std::endl;
	std::cerr << "Read ROOT file" << std::endl;
	std::cerr << "Visualizer for 0.1c NIuTPC" << std::endl;
	std::cerr << "Version 0.1" << std::endl;
	std::cerr << "======================================" << std::endl;
	std::cerr << "input file name   : " << dirfilename << endl;
	std::cerr << "input config file : " << conffilename << endl;
	ni_conf->PrintConfigJSON();
	std::cerr << "======================================" << std::endl;
	TApplication app("app",&argc,argv);  

	TFile *f=new TFile(dirfilename.c_str());
	if(!f){
		std::cerr << "ERROR: Cant find file" << std::endl;
		return 1;
	}
	TTree *tree=(TTree*)f->Get("Tree");

	Event *event = new Event();
	int nevent=tree->GetEntries();
	std::cerr << "NumberOfEvents : " << nevent << std::endl;
	
	tree->SetBranchAddress("Event",&event);
	
	// Get Event Info
	tree->GetEntry(ev_num);
	std::cerr << "Draw Event Number : " << ev_num << std::endl;
	int trigger_num = event->GetTrigger();
	
	vector<vector<double>> a_hg_adc = event->GetAnodeHGADC();
	vector<vector<double>> a_lg_adc = event->GetAnodeLGADC();
	vector<vector<double>> c_hg_adc = event->GetCathodeHGADC();
	vector<vector<double>> c_lg_adc = event->GetCathodeLGADC();
	
	// Pedestal calc
	double pedestal_a_hg[64];
	double pedestal_a_lg[64];
	double pedestal_c_hg[64];
	double pedestal_c_lg[64];
	double noise_a_hg[64];
	double noise_a_lg[64];
	double noise_c_hg[64];
	double noise_c_lg[64];
	bool   mask_a_hg[64];
	bool   mask_a_lg[64];
	bool   mask_c_hg[64];
	bool   mask_c_lg[64];
	for(int j=0;j<N_CHANNEL;j++){
		pedestal_a_hg[j]=0;
		pedestal_a_lg[j]=0;
		pedestal_c_hg[j]=0;
		pedestal_c_lg[j]=0;
		noise_a_hg[j]=0;
		noise_a_lg[j]=0;
		noise_c_hg[j]=0;
		noise_c_lg[j]=0;
		for(int k=0;k<offset_sampling;k++){
			pedestal_a_hg[j] += a_hg_adc[j][k]/offset_sampling;
			pedestal_a_lg[j] += a_lg_adc[j][k]/offset_sampling;
			pedestal_c_hg[j] += c_hg_adc[j][k]/offset_sampling;
			pedestal_c_lg[j] += c_lg_adc[j][k]/offset_sampling;
		}
		for(int k=0;k<offset_sampling;k++){
            noise_a_hg[j]+=pow(a_hg_adc[j][k]-pedestal_a_hg[j], 2)/offset_sampling;
            noise_a_lg[j]+=pow(a_lg_adc[j][k]-pedestal_a_lg[j], 2)/offset_sampling;
            noise_c_hg[j]+=pow(c_hg_adc[j][k]-pedestal_c_hg[j], 2)/offset_sampling;
            noise_c_lg[j]+=pow(c_lg_adc[j][k]-pedestal_c_lg[j], 2)/offset_sampling;
		}
        noise_a_hg[j]=sqrt(noise_a_hg[j]);
        noise_a_lg[j]=sqrt(noise_a_lg[j]);
        noise_c_hg[j]=sqrt(noise_c_hg[j]);
        noise_c_lg[j]=sqrt(noise_c_lg[j]);

        mask_a_hg[j] = ( noise_a_hg[j] > 50.0 ) ? true : false;
        mask_a_lg[j] = ( noise_a_lg[j] > 5.0   ) ? true : false;
        mask_c_hg[j] = ( noise_c_hg[j] > 50.0 ) ? true : false;
        mask_c_lg[j] = ( noise_c_lg[j] > 5.0   ) ? true : false;

        // if( mask_a_hg[j] == true ) std::cout << "noise_a_hg[" << j << "]: " << noise_a_hg[j] << std::endl; 
        // if( mask_a_lg[j] == true ) std::cout << "noise_a_lg[" << j << "]: " << noise_a_lg[j] << std::endl; 
        // if( mask_c_hg[j] == true ) std::cout << "noise_c_hg[" << j << "]: " << noise_c_hg[j] << std::endl; 
        // if( mask_c_lg[j] == true ) std::cout << "noise_c_lg[" << j << "]: " << noise_c_lg[j] << std::endl; 
	}


	//++++++++++++++++++++++++++++++++++++++++++
	//  Fill waveform
	//++++++++++++++++++++++++++++++++++++++++++
	//Even CH -> High Gain
	//Odd CH  -> Low  Gain
	TGraph* g_a_hg[N_CHANNEL];
	TGraph* g_a_lg[N_CHANNEL];
	TGraph* g_c_hg[N_CHANNEL];
	TGraph* g_c_lg[N_CHANNEL];
	for(int i=0;i<N_CHANNEL;i++){
		g_a_hg[i] = new TGraph();
		g_a_lg[i] = new TGraph();
		g_c_hg[i] = new TGraph();
		g_c_lg[i] = new TGraph();
	}


	TH2F *h_strip_a_hg=new TH2F("h_strip_a_hg","h_strip_a_hg",SAMPLING_NUM,0,SAMPLING_NUM,N_CHANNEL,0,N_CHANNEL);
	TH2F *h_strip_a_lg=new TH2F("h_strip_a_lg","h_strip_a_lg",SAMPLING_NUM,0,SAMPLING_NUM,N_CHANNEL,0,N_CHANNEL);
	TH2F *h_strip_c_hg=new TH2F("h_strip_c_hg","h_strip_c_hg",SAMPLING_NUM,0,SAMPLING_NUM,N_CHANNEL,0,N_CHANNEL);
	TH2F *h_strip_c_lg=new TH2F("h_strip_c_lg","h_strip_c_lg",SAMPLING_NUM,0,SAMPLING_NUM,N_CHANNEL,0,N_CHANNEL);
	h_strip_a_hg->SetStats(0);
	h_strip_a_lg->SetStats(0);
	h_strip_c_hg->SetStats(0);
	h_strip_c_lg->SetStats(0);
	for(int j=0;j<N_CHANNEL;j++){
		for(int k=0;k<SAMPLING_NUM;k++){
			g_a_hg[j]->SetPoint(k, k/SAMPLING_HELZ*1e6, a_hg_adc.at(j).at(k)-pedestal_a_hg[j] - (63-j)*100);
			g_a_lg[j]->SetPoint(k, k/SAMPLING_HELZ*1e6, a_lg_adc.at(j).at(k)-pedestal_a_lg[j] - (63-j)*10);
			g_c_hg[j]->SetPoint(k, k/SAMPLING_HELZ*1e6, c_hg_adc.at(j).at(k)-pedestal_c_hg[j] - (63-j)*100);
			g_c_lg[j]->SetPoint(k, k/SAMPLING_HELZ*1e6, c_lg_adc.at(j).at(k)-pedestal_c_lg[j] - (63-j)*10);
            g_a_hg[j]->SetMarkerColor(MyPalette[j]);
            g_a_lg[j]->SetMarkerColor(MyPalette[j]);
            g_c_hg[j]->SetMarkerColor(MyPalette[j]);
            g_c_lg[j]->SetMarkerColor(MyPalette[j]);

            double a_hg_adc_abs = mask_a_hg[j] ? 0.0 : a_hg_adc[j][k]-pedestal_a_hg[j];
            double a_lg_adc_abs = mask_a_lg[j] ? 0.0 : a_lg_adc[j][k]-pedestal_a_lg[j];
            double c_hg_adc_abs = mask_c_hg[j] ? 0.0 : c_hg_adc[j][k]-pedestal_c_hg[j];
            double c_lg_adc_abs = mask_c_lg[j] ? 0.0 : c_lg_adc[j][k]-pedestal_c_lg[j];
			h_strip_a_hg->SetBinContent(k+1,j+1,a_hg_adc_abs);
			h_strip_a_lg->SetBinContent(k+1,j+1,a_lg_adc_abs);
			h_strip_c_hg->SetBinContent(k+1,j+1,-c_hg_adc_abs);
			h_strip_c_lg->SetBinContent(k+1,j+1,-c_lg_adc_abs);
		}
	}
	
	//++++++++++++++++++++++++++++++++++++++++++
	//  Main Peak analysis
	//++++++++++++++++++++++++++++++++++++++++++
	//anode
	TGraph* g_a_hg_main_peak = new TGraph();
	TGraph* g_a_lg_main_peak = new TGraph();
	TGraph* g_a_main_rise = new TGraph();
	TGraph* g_a_main_fall = new TGraph();
	TH1D* h_a_hg_charge = new TH1D("h_a_hg_charge","h_a_hg_charge",N_CHANNEL,0,N_CHANNEL);
	TH1D* h_a_lg_charge = new TH1D("h_a_lg_charge","h_a_lg_charge",N_CHANNEL,0,N_CHANNEL);
	double sum_a_hg_pulse_height = 0;
	double sum_a_lg_pulse_height = 0;
	double sum_a_charge = 0;
	vector<double> hg_a_pulse_height;
	vector<double> lg_a_pulse_height;
	vector<double> hg_a_charge;
	vector<double> lg_a_charge;
	vector<double> hg_a_peak_time;
	vector<double> lg_a_peak_time;
	//cathode
	TGraph* g_c_hg_main_peak = new TGraph();
	TGraph* g_c_lg_main_peak = new TGraph();
	TGraph* g_c_main_rise = new TGraph();
	TGraph* g_c_main_fall = new TGraph();
	TH1D* h_c_hg_charge = new TH1D("h_c_hg_charge","h_c_hg_charge",N_CHANNEL,0,N_CHANNEL);
	TH1D* h_c_lg_charge = new TH1D("h_c_lg_charge","h_c_lg_charge",N_CHANNEL,0,N_CHANNEL);
	double sum_c_hg_pulse_height = 0;
	double sum_c_lg_pulse_height = 0;
	double sum_c_charge = 0;
	vector<double> hg_c_pulse_height;
	vector<double> lg_c_pulse_height;
	vector<double> hg_c_charge;
	vector<double> lg_c_charge;
	vector<double> hg_c_peak_time;
	vector<double> lg_c_peak_time;;
	//minority
	TGraph* g_mino_peak = new TGraph();
	TGraph* g_mino_search = new TGraph();
	int mino_search_point = 0;

	for(int j=0;j<N_CHANNEL;j++){
		//peak saerch
		//anode
		double hg_a_mainpeak_time=-1;
		double lg_a_mainpeak_time=-1;
		double hg_a_mainrise_time=-1;
		double lg_a_mainrise_time=-1;
		double hg_a_mainfall_time=-1;
		double lg_a_mainfall_time=-1;
		double hg_a_pulse_max=-9999;
		double lg_a_pulse_max=-9999;
		double hg_a_this_charge = 0;
		double lg_a_this_charge = 0;
		int hg_a_rise_flag=0;
		//cathode
		double hg_c_mainpeak_time=-1;
		double lg_c_mainpeak_time=-1;
		double hg_c_mainrise_time=-1;
		double lg_c_mainrise_time=-1;
		double hg_c_mainfall_time=-1;
		double lg_c_mainfall_time=-1;
		double hg_c_pulse_max=9999;
		double lg_c_pulse_max=9999;
		double hg_c_this_charge = 0;
		double lg_c_this_charge = 0;
		int hg_c_rise_flag=0;
		//minority
		TH1D* h_mino_search = new TH1D("h_mino_search","h_mino_search",SAMPLING_NUM,0,SAMPLING_NUM/SAMPLING_HELZ*1e6);
		for(int k=0;k<SAMPLING_NUM;k++){
			//anode
			if(a_hg_adc.at(j).at(k)-pedestal_a_hg[j] > tot_anode_threshold){
				if(hg_a_rise_flag==0){
					hg_a_rise_flag=1;
					hg_a_mainrise_time = k/SAMPLING_HELZ*1e6;
				}
				hg_a_mainfall_time = k/SAMPLING_HELZ*1e6;
				hg_a_this_charge += a_hg_adc.at(j).at(k)-pedestal_a_hg[j];
				lg_a_this_charge += a_lg_adc.at(j).at(k)-pedestal_a_lg[j];
			}
			if(hg_a_pulse_max < a_hg_adc.at(j).at(k)-pedestal_a_hg[j]){ //use HG
				hg_a_pulse_max = a_hg_adc.at(j).at(k)-pedestal_a_hg[j];
				if(hg_a_pulse_max>tot_anode_threshold){
					hg_a_mainpeak_time=k/SAMPLING_HELZ*1e6;
				}
			}
			if(lg_a_pulse_max < a_lg_adc.at(j).at(k)-pedestal_a_lg[j]){ //use LG
				lg_a_pulse_max = a_lg_adc.at(j).at(k)-pedestal_a_lg[j];
				if(lg_a_pulse_max>lg_anode_threshold){
					lg_a_mainpeak_time=k/SAMPLING_HELZ*1e6;
				}
			}
			//cathode
			if(c_hg_adc.at(j).at(k)-pedestal_c_hg[j] < tot_cathode_threshold){
				if(hg_c_rise_flag==0){
					hg_c_rise_flag=1;
					hg_c_mainrise_time = k/SAMPLING_HELZ*1e6;
				}
				hg_c_mainfall_time = k/SAMPLING_HELZ*1e6;
				hg_c_this_charge += c_hg_adc.at(j).at(k)-pedestal_c_hg[j];
				lg_c_this_charge += c_lg_adc.at(j).at(k)-pedestal_c_lg[j];
			}
			if(hg_c_pulse_max > c_hg_adc.at(j).at(k)-pedestal_c_hg[j]){ //use HG
				hg_c_pulse_max = c_hg_adc.at(j).at(k)-pedestal_c_hg[j];
				if(hg_c_pulse_max<tot_cathode_threshold){
					hg_c_mainpeak_time=k/SAMPLING_HELZ*1e6;
				}
			}
			if(lg_c_pulse_max > c_lg_adc.at(j).at(k)-pedestal_c_lg[j]){ //use LG
				lg_c_pulse_max = c_lg_adc.at(j).at(k)-pedestal_c_lg[j];
				if(lg_c_pulse_max<lg_cathode_threshold){
					lg_c_mainpeak_time=k/SAMPLING_HELZ*1e6;
				}
			}
			h_mino_search->SetBinContent(k,a_hg_adc[j][k]-pedestal_a_hg[j]);//anode only
		}
		//++++++++++++++++++++++++++++++++++++++++++
		//  Minority Peak analysis
		//++++++++++++++++++++++++++++++++++++++++++
		double this_mino_time=-1;
		int ROI_min_bin = h_mino_search->FindBin(lg_a_mainpeak_time - minority_ROI_end);//us
		int ROI_max_bin = h_mino_search->FindBin(lg_a_mainpeak_time - minority_ROI_start);//us
		int ROI_bins = ROI_max_bin - ROI_min_bin + 1;
		if(lg_a_mainpeak_time==-1){
			//skip
		}else{
			TH1D* h_ROI = new TH1D("h_ROI","h_ROI",ROI_bins,lg_a_mainpeak_time-minority_ROI_end,lg_a_mainpeak_time-minority_ROI_start);
			for(int k=0;k<ROI_max_bin;k++){
				h_ROI->SetBinContent(k,h_mino_search->GetBinContent(ROI_min_bin+k));
			}
			TSpectrum* s = new TSpectrum(1000);
			int find_mino_peak = s->Search(h_ROI,1,"nodraw",0.5);
			double* some_mino_time   = s->GetPositionX();
			double* some_mino_height = s->GetPositionY();
			double mino_peak_max =0;
			double mino_dt_min = minority_ROI_end;
			for(int pks=0;pks<find_mino_peak;pks++){
				if(minority_threshold > some_mino_height[pks])continue;
				//maxmum method
				/*
                  if(some_mino_height[pks]>mino_peak_max){
                  g_mino_search->SetPoint(mino_search_point,some_mino_time[pks],some_mino_height[pks]-(63-j)*50);
                  mino_search_point++;
                  this_mino_time = some_mino_time[pks];
                  mino_peak_max = some_mino_height[pks];
                  }
				*/
				if(lg_a_mainpeak_time-some_mino_time[pks]<mino_dt_min){
					g_mino_search->SetPoint(mino_search_point,some_mino_time[pks],some_mino_height[pks]-(63-j)*50);
					mino_search_point++;
					this_mino_time = some_mino_time[pks];
					mino_dt_min = lg_a_mainpeak_time-some_mino_time[pks]; 
				}
			}
			s->Delete();
			h_ROI->Delete();
		}
		h_mino_search->Delete();
		//anode
		g_a_hg_main_peak->SetPoint(j,hg_a_mainpeak_time*10/4,j+0.5);
		g_a_main_rise->SetPoint(j,hg_a_mainrise_time*10/4,j+0.5);
		g_a_main_fall->SetPoint(j,hg_a_mainfall_time*10/4,j+0.5);
		g_a_lg_main_peak->SetPoint(j,lg_a_mainpeak_time*10/4,j+0.5);
		h_a_hg_charge->SetBinContent(j,hg_a_this_charge);
		h_a_lg_charge->SetBinContent(j,lg_a_this_charge);
		//cathode
		g_c_hg_main_peak->SetPoint(j,hg_c_mainpeak_time*10/4,j+0.5);
		g_c_main_rise->SetPoint(j,hg_c_mainrise_time*10/4,j+0.5);
		g_c_main_fall->SetPoint(j,hg_c_mainfall_time*10/4,j+0.5);
		g_c_lg_main_peak->SetPoint(j,lg_c_mainpeak_time*10/4,j+0.5);
		h_c_hg_charge->SetBinContent(j,hg_c_this_charge);
		h_c_lg_charge->SetBinContent(j,lg_c_this_charge);
		//minority
		g_mino_peak->SetPoint(j,this_mino_time*10/4,j+0.5);
	}



	//++++++++++++++++++++++++++++++++++++++++++
	//  Draw
	//++++++++++++++++++++++++++++++++++++++++++
	gStyle->SetPalette(kRainBow);
	TCanvas *c_strip=new TCanvas("c_strip","",800,800);
	c_strip->Divide(2,2);
	
	c_strip->cd(1);
	g_a_hg_main_peak->SetMarkerStyle(8);
	g_a_hg_main_peak->SetMarkerSize(0.3);
	g_a_hg_main_peak->SetMarkerColor(kMagenta);
	g_a_main_rise->SetMarkerStyle(8);
	g_a_main_rise->SetMarkerSize(0.3);
	g_a_main_rise->SetMarkerColor(kRed);
	g_a_main_fall->SetMarkerStyle(8);
	g_a_main_fall->SetMarkerSize(0.3);
	g_a_main_fall->SetMarkerColor(kGreen);
	g_mino_peak->SetMarkerStyle(8);
	g_mino_peak->SetMarkerSize(0.3);
	g_mino_peak->SetMarkerColor(kPink);
	h_strip_a_hg->GetXaxis()->SetTitle("clock(1clock=0.4us)");
	h_strip_a_hg->GetYaxis()->SetTitle("CH");
	h_strip_a_hg->Draw("colz");
	// g_a_hg_main_peak->Draw("same p");
	// g_a_main_rise->Draw("same p");
	// g_a_main_fall->Draw("same p");
	// g_mino_peak->Draw("same p");
	
	c_strip->cd(2);
	g_a_lg_main_peak->SetMarkerStyle(8);
	g_a_lg_main_peak->SetMarkerSize(0.3);
	g_a_lg_main_peak->SetMarkerColor(kMagenta);
	h_strip_a_lg->GetXaxis()->SetTitle("clock(1clock=0.4us)");
	h_strip_a_lg->GetYaxis()->SetTitle("CH");
	h_strip_a_lg->GetZaxis()->SetRangeUser(-10.0, 10.0);
	h_strip_a_lg->Draw("colz");
	// g_a_lg_main_peak->Draw("same p");
	
	c_strip->cd(3);
	g_c_hg_main_peak->SetMarkerStyle(8);
	g_c_hg_main_peak->SetMarkerSize(0.3);
	g_c_hg_main_peak->SetMarkerColor(kMagenta);
	g_c_main_rise->SetMarkerStyle(8);
	g_c_main_rise->SetMarkerSize(0.3);
	g_c_main_rise->SetMarkerColor(kRed);
	g_c_main_fall->SetMarkerStyle(8);
	g_c_main_fall->SetMarkerSize(0.3);
	g_c_main_fall->SetMarkerColor(kGreen);
	h_strip_c_hg->GetXaxis()->SetTitle("clock(1clock=0.4us)");
	h_strip_c_hg->GetYaxis()->SetTitle("CH");
	h_strip_c_hg->Draw("colz");
	// g_c_hg_main_peak->Draw("same p");
	// g_c_main_rise->Draw("same p");
	// g_c_main_fall->Draw("same p");
	
	c_strip->cd(4);
	g_c_lg_main_peak->SetMarkerStyle(8);
	g_c_lg_main_peak->SetMarkerSize(0.3);
	g_c_lg_main_peak->SetMarkerColor(kMagenta);
	h_strip_c_lg->GetXaxis()->SetTitle("clock(1clock=0.4us)");
	h_strip_c_lg->GetYaxis()->SetTitle("CH");
	h_strip_c_lg->Draw("colz");
	// g_c_lg_main_peak->Draw("same p");

    c_strip->SaveAs( "recoil.png" );
    // c_strip->SaveAs( Form( "output6/strip_%d.png", ev_num ) );
	
	TCanvas *c_wave=new TCanvas("c_wave","",1000,1000);
	c_wave->Divide(2,2);
	c_wave->cd(1)->DrawFrame(0,-4000,1600,1000,"HG waveform;us(4000sampling);mV");
	c_wave->cd(2)->DrawFrame(0,-4000,1600,1000,"HG waveform;us(4000sampling);mV");
	c_wave->cd(3)->DrawFrame(0,-400,1600,100,"LG waveform;us(4000sampling);mV");
	c_wave->cd(4)->DrawFrame(0,-400,1600,100,"LG waveform;us(4000sampling);mV");
	for(int i=0;i<N_CHANNEL;i++){
		c_wave->cd(1); if( mask_a_hg[i] == false ) g_a_hg[i]->Draw("p same");
		c_wave->cd(3); if( mask_a_lg[i] == false ) g_a_lg[i]->Draw("p same");
		c_wave->cd(2); if( mask_c_hg[i] == false ) g_c_hg[i]->Draw("p same");
		c_wave->cd(4); if( mask_c_lg[i] == false ) g_c_lg[i]->Draw("p same");
	}
	c_wave->cd(1);
	g_mino_search->SetMarkerStyle(8);
	g_mino_search->SetMarkerSize(0.5);
	g_mino_search->SetMarkerColor(kGreen+3);
	// g_mino_search->Draw("same p");

	TCanvas *c_charge = new TCanvas("c_charge","c_charge",1000,1000);
	c_charge->Divide(2,2);

	c_charge->cd(1);
	h_a_hg_charge->SetFillColor(kAzure+10);
	h_a_hg_charge->SetLineColor(kAzure+1);
	h_a_hg_charge->SetLineWidth(2);
	h_a_hg_charge->GetXaxis()->SetTitle("CH");	
	h_a_hg_charge->GetYaxis()->SetTitle("charge(mV*0.4us)");	
	h_a_hg_charge->Draw("hist");	
	c_charge->cd(2);
	h_c_hg_charge->SetFillColor(kMagenta-10);
	h_c_hg_charge->SetLineColor(kMagenta+1);
	h_c_hg_charge->SetLineWidth(2);
	h_c_hg_charge->GetXaxis()->SetTitle("CH");	
	h_c_hg_charge->GetYaxis()->SetTitle("charge(mV*0.4us)");	
	h_c_hg_charge->Draw("hist");	

	c_charge->cd(3);
	h_a_lg_charge->SetFillColor(kAzure+10);
	h_a_lg_charge->SetLineColor(kAzure+1);
	h_a_lg_charge->SetLineWidth(2);
	h_a_lg_charge->GetXaxis()->SetTitle("CH");	
	h_a_lg_charge->GetYaxis()->SetTitle("charge(mV*0.4us)");	
	h_a_lg_charge->Draw("hist");	
	c_charge->cd(4);
	h_c_lg_charge->SetFillColor(kMagenta-10);
	h_c_lg_charge->SetLineColor(kMagenta+1);
	h_c_lg_charge->SetLineWidth(2);
	h_c_lg_charge->GetXaxis()->SetTitle("CH");	
	h_c_lg_charge->GetYaxis()->SetTitle("charge(mV*0.4us)");	
	h_c_lg_charge->Draw("hist");	

	//++++++++++++++++++++++++++++++++++++++++++
	//  Output
	//++++++++++++++++++++++++++++++++++++++++++
    /*
	std::string filename;
	std::string::size_type pos = dirfilename.find("/");
	while(pos != std::string::npos){
		dirfilename = dirfilename.substr(pos+1);
		pos = dirfilename.find("/");
	}
	pos = dirfilename.find(".root");
	filename = dirfilename.substr(0,pos);
	std::ostringstream oss_ev_num;
	oss_ev_num << ev_num;
	std::string outfilename = filename.substr(0, index) + "_event" + oss_ev_num.str() + ".root";
	TFile *fout=new TFile(outfilename.c_str(),"recreate");

    // TTree* event_tree = new TTree("event_tree","event_tree");
    // int outtree_ev_num;
    // vector<double> outtree_hg_a_pulse_height;
    // event_tree->Branch("ev_num",ev_num);

	for(int i=0;i<N_CHANNEL;i++){
		g_a_hg[i]->Write();
		g_a_lg[i]->Write();
		g_c_hg[i]->Write();
		g_c_lg[i]->Write();
	}
	h_strip_a_hg->Write();
	h_strip_a_lg->Write();
	h_strip_c_hg->Write();
	h_strip_c_lg->Write();
	fout->Close();
    */
	app.Run();

	return 1;
}
