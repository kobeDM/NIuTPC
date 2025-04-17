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
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TClassTable.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TSpectrum.h"
#include "TLatex.h"
//USER
#include "Event.h"
#include "NIConfig.hh"
#include "NAPStyle.h"

#define N_CHANNEL 64
#define N_STRIP 32
#define SAMPLING_NUM 4000
#define SAMPLING_HELZ 2.5e6 //Hz

#define DEBUG 0
using namespace std;

TLatex* CreateDrawText( const double&       x,
                        const double&       y,
                        const std::string&  text,
                        const double&       size = 0.05,
                        const Color_t&      color = 1 )
{
    if( text.size( ) <= 0 ) return nullptr;
    TLatex l;
    l.SetNDC( );
    l.SetTextColor( color );
    l.SetTextSize( size );
    return l.DrawLatex( x, y, text.c_str( ) );
}


//------------- main -------------------------------------------------
int main(int argc,char *argv[]){

    SetNAPStyle( );

	Int_t MyPalette[64];
	Double_t Red   [] = {1.0, 0.0, 1.0};
	Double_t Green [] = {0.0, 0.0, 0.0};
	Double_t Blue  [] = {0.0, 1.0, 1.0};
	Double_t Length[] = {0.0, 0.5, 1.0};
	Int_t FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 64);
	for (int i=0;i<64;i++) MyPalette[i] = FI+i;

	if(argc <4){
		std::cerr << "Usage:" << std::endl;
		std::cerr << "/.NIEventShow [filename.root] [config.json] [event number]" << std::endl;
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
	double hg_anode_threshold    = ni_conf->hg_anode_threshold;
	double lg_anode_threshold    = ni_conf->lg_anode_threshold;
	double hg_cathode_threshold  = ni_conf->hg_cathode_threshold;
	double lg_cathode_threshold  = ni_conf->lg_cathode_threshold;
	double minority_threshold    = ni_conf->minority_threshold;
	double minority_ROI_range    = ni_conf->minority_ROI_range;
	double minority_ROI_offset   = ni_conf->minority_ROI_offset;

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
	// TApplication app("app",&argc,argv);  

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

        // mask_a_hg[j] = ( noise_a_hg[j] > 50.0 ) ? true : false;
        // mask_a_lg[j] = ( noise_a_lg[j] > 5.0  ) ? true : false;
        // mask_c_hg[j] = ( noise_c_hg[j] > 50.0 ) ? true : false;
        // mask_c_lg[j] = ( noise_c_lg[j] > 5.0  ) ? true : false;


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

    TGraph* g_a_hg_ave = new TGraph();
    TGraph* g_c_hg_ave = new TGraph();
    TGraph* g_a_lg_ave = new TGraph();
    TGraph* g_c_lg_ave = new TGraph();

	TH2F *h_strip_a_hg=new TH2F("h_strip_a_hg","h_strip_a_hg",SAMPLING_NUM,0,SAMPLING_NUM,N_CHANNEL,0,N_CHANNEL);
	TH2F *h_strip_a_lg=new TH2F("h_strip_a_lg","h_strip_a_lg",SAMPLING_NUM,0,SAMPLING_NUM,N_CHANNEL,0,N_CHANNEL);
	TH2F *h_strip_c_hg=new TH2F("h_strip_c_hg","h_strip_c_hg",SAMPLING_NUM,0,SAMPLING_NUM,N_CHANNEL,0,N_CHANNEL);
	TH2F *h_strip_c_lg=new TH2F("h_strip_c_lg","h_strip_c_lg",SAMPLING_NUM,0,SAMPLING_NUM,N_CHANNEL,0,N_CHANNEL);
	h_strip_a_hg->SetStats(0);
	h_strip_a_lg->SetStats(0);
	h_strip_c_hg->SetStats(0);
	h_strip_c_lg->SetStats(0);

	// TH2F *h_track_a_lg=new TH2F("h_track_a_lg","h_track_a_lg",SAMPLING_NUM,-1000,1000,N_CHANNEL,-10,10);
	// TH2F *h_track_c_lg=new TH2F("h_track_c_lg","h_track_c_lg",SAMPLING_NUM,0,SAMPLING_NUM/SAMPLING_HELZ*1e3*driftV_main,N_CHANNEL,0,N_CHANNEL*0.04);
	TH2F *h_track_a_lg=new TH2F("h_track_a_lg","h_track_a_lg",SAMPLING_NUM,0,SAMPLING_NUM/SAMPLING_HELZ*1e6,N_CHANNEL,N_CHANNEL*0.04*-0.5,N_CHANNEL*0.04*0.5);
	TH2F *h_track_c_lg=new TH2F("h_track_c_lg","h_track_c_lg",SAMPLING_NUM,0,SAMPLING_NUM/SAMPLING_HELZ*1e6,N_CHANNEL,N_CHANNEL*0.04*-0.5,N_CHANNEL*0.04*0.5);

    for(int k=0;k<SAMPLING_NUM;k++){
        
        double a_hg_sum = 0.0, c_hg_sum = 0.0;
        double a_lg_sum = 0.0, c_lg_sum = 0.0;
        int averaged_ch = 0;
        for(int j=0;j<N_CHANNEL;j++){
			g_a_hg[j]->SetPoint(k, k/SAMPLING_HELZ*1e6, a_hg_adc.at(j).at(k)-pedestal_a_hg[j] - (31-j)*100);
			g_a_lg[j]->SetPoint(k, k/SAMPLING_HELZ*1e6, a_lg_adc.at(j).at(k)-pedestal_a_lg[j] - (31-j)*10);
			g_c_hg[j]->SetPoint(k, k/SAMPLING_HELZ*1e6, c_hg_adc.at(j).at(k)-pedestal_c_hg[j] - (31-j)*100);
			g_c_lg[j]->SetPoint(k, k/SAMPLING_HELZ*1e6, c_lg_adc.at(j).at(k)-pedestal_c_lg[j] - (31-j)*10);
            g_a_hg[j]->SetMarkerColor(MyPalette[j]);
            g_a_lg[j]->SetMarkerColor(MyPalette[j]);
            g_c_hg[j]->SetMarkerColor(MyPalette[j]);
            g_c_lg[j]->SetMarkerColor(MyPalette[j]);

            g_a_hg[j]->SetMarkerStyle( 1 );
            g_a_lg[j]->SetMarkerStyle( 1 );
            g_c_hg[j]->SetMarkerStyle( 1 );
            g_c_lg[j]->SetMarkerStyle( 1 );

            double a_hg_adc_abs = mask_a_hg[j] ? 0.0 : a_hg_adc[j][k]-pedestal_a_hg[j];
            double a_lg_adc_abs = mask_a_lg[j] ? 0.0 : a_lg_adc[j][k]-pedestal_a_lg[j];
            double c_hg_adc_abs = mask_c_hg[j] ? 0.0 : c_hg_adc[j][k]-pedestal_c_hg[j];
            double c_lg_adc_abs = mask_c_lg[j] ? 0.0 : c_lg_adc[j][k]-pedestal_c_lg[j];
			h_strip_a_hg->SetBinContent(k+1,j+1,a_hg_adc_abs);
			h_strip_a_lg->SetBinContent(k+1,j+1,a_lg_adc_abs);
			h_strip_c_hg->SetBinContent(k+1,j+1,-c_hg_adc_abs);
			h_strip_c_lg->SetBinContent(k+1,j+1,-c_lg_adc_abs);

			h_track_a_lg->SetBinContent( k+1 ,N_CHANNEL-j,a_lg_adc_abs);
			h_track_c_lg->SetBinContent( k+1 ,N_CHANNEL-j,-c_lg_adc_abs);
			// h_track_a_lg->SetBinContent( k+1 ,j+1,a_lg_adc_abs);
			// h_track_c_lg->SetBinContent( k+1 ,j+1,-c_lg_adc_abs);

            if( mask_a_hg[j] == false ) {
                a_hg_sum += a_hg_adc.at(j).at(k)-pedestal_a_hg[j];
                a_lg_sum += a_lg_adc.at(j).at(k)-pedestal_a_lg[j];

                c_hg_sum += c_hg_adc.at(j).at(k)-pedestal_c_hg[j];
                c_lg_sum += c_lg_adc.at(j).at(k)-pedestal_c_lg[j];
                averaged_ch += 1;
            }
		}

        // g_a_hg_ave->SetPoint( k, k/SAMPLING_HELZ*1e6, a_hg_sum / ((double)N_CHANNEL ) );
        // g_c_hg_ave->SetPoint( k, k/SAMPLING_HELZ*1e6, c_hg_sum / ((double)N_CHANNEL ) );
        // g_a_lg_ave->SetPoint( k, k/SAMPLING_HELZ*1e6, a_lg_sum / ((double)N_CHANNEL ) );
        // g_c_lg_ave->SetPoint( k, k/SAMPLING_HELZ*1e6, c_lg_sum / ((double)N_CHANNEL ) );

        g_a_hg_ave->SetPoint( k, k/SAMPLING_HELZ*1e6, a_hg_sum / ((double)averaged_ch ) );
        g_c_hg_ave->SetPoint( k, k/SAMPLING_HELZ*1e6, c_hg_sum / ((double)averaged_ch ) );
        g_a_lg_ave->SetPoint( k, k/SAMPLING_HELZ*1e6, a_lg_sum / ((double)averaged_ch ) );
        g_c_lg_ave->SetPoint( k, k/SAMPLING_HELZ*1e6, c_lg_sum / ((double)averaged_ch ) );

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
	TGraph* g_mino_search_max = new TGraph();
	int mino_search_point = 0;

	TGraph* g_main_search = new TGraph();
	int main_search_point = 0;

	TH1D* h_dt = new TH1D("h_dt","h_dt", 100, 0.0, 200.0 );

    std::map< int, TGraph* > hg_minoAnaWfTable;
    std::map< int, TGraph* > lg_minoAnaWfTable;

    TH1D* h_all_mino_wf_hg = new TH1D("h_all_mino_wf_hg","h_all_mino_wf_hg",1024, -204.8, 204.8 );
    TH1D* h_all_mino_wf_lg = new TH1D("h_all_mino_wf_lg","h_all_mino_wf_lg",1024, -204.8, 204.8 );

	for(int j=0;j<N_CHANNEL;j++){
		//peak saerch
		//anode
		double hg_a_mainpeak_time=-1;
		double lg_a_mainpeak_time=-1;
        int    hg_a_mainpeak_idx =-1;
        int    lg_a_mainpeak_idx =-1;
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
			if(a_hg_adc.at(j).at(k)-pedestal_a_hg[j] > hg_anode_threshold){
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
				if(hg_a_pulse_max>hg_anode_threshold){
					hg_a_mainpeak_time=k/SAMPLING_HELZ*1e6;
                    hg_a_mainpeak_idx=k;
				}
			}
			if(lg_a_pulse_max < a_lg_adc.at(j).at(k)-pedestal_a_lg[j]){ //use LG
				lg_a_pulse_max = a_lg_adc.at(j).at(k)-pedestal_a_lg[j];
				if(lg_a_pulse_max>lg_anode_threshold){
					lg_a_mainpeak_time=k/SAMPLING_HELZ*1e6;
                    lg_a_mainpeak_idx=k;
				}
			}
			//cathode
			if(c_hg_adc.at(j).at(k)-pedestal_c_hg[j] < hg_cathode_threshold){
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
				if(hg_c_pulse_max<hg_cathode_threshold){
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
        double ROI_min = lg_a_mainpeak_time - minority_ROI_offset - minority_ROI_range;
        double ROI_max = lg_a_mainpeak_time - minority_ROI_offset;
		int ROI_min_bin = h_mino_search->FindBin(ROI_min);//us
		int ROI_max_bin = h_mino_search->FindBin(ROI_max);//us
		int ROI_bins = ROI_max_bin - ROI_min_bin + 1;
        std::cout << "ROI min: " << ROI_min << "\t" << ROI_max << "\t" << lg_a_mainpeak_time << std::endl;


		if(lg_a_mainpeak_time==-1){
			//skip
		}else{
			TH1D* h_ROI = new TH1D("h_ROI","h_ROI",ROI_bins,ROI_min,ROI_max);
			for(int k=0;k<ROI_max_bin;k++){
				h_ROI->SetBinContent(k,h_mino_search->GetBinContent(ROI_min_bin+k));
			}
			TSpectrum* s = new TSpectrum(1000);
			int find_mino_peak = s->Search(h_ROI,1,"nodraw",0.5);
			double* some_mino_time   = s->GetPositionX();
			double* some_mino_height = s->GetPositionY();
			double mino_peak_max =-1.0;
			double mino_dt_min = ROI_min;
			// for(int pks=0;pks<find_mino_peak;pks++){
            //      if( j == 0 || j == 1 || j == 2 ) continue; // need to fix
			// 	if(minority_threshold > some_mino_height[pks])continue;
				
            //     g_mino_search->SetPoint(mino_search_point,some_mino_time[pks],some_mino_height[pks]-(31-j)*100);
            //     if(some_mino_height[pks]>mino_peak_max){
            //         g_mino_search_max->SetPoint(mino_search_point,some_mino_time[pks],some_mino_height[pks]-(31-j)*100);
            //         std::cout << j << "\t" << some_mino_time[pks] << "\t" << some_mino_height[pks] << std::endl;
            //         mino_search_point++;
            //         this_mino_time = some_mino_time[pks];
            //         mino_peak_max = some_mino_height[pks];
            //     }

            //     double eachDt = lg_a_mainpeak_time - some_mino_time[pks];
            //     h_dt->Fill( eachDt );
			// }

            std::cout << "test000" << std::endl;
            // if( mino_peak_max > 0.0 ) {
            //     TGraph* pHgGraph = new TGraph( );
            //     TGraph* pLgGraph = new TGraph( );
            //     TH1D histHG("histHG","histHG",1024, -204.8, 204.8);
            //     TH1D histLG("histLG","histLG",1024, -204.8, 204.8);
            //     for( int relIdx = 0; relIdx < 1024; ++relIdx ) {
            //         pHgGraph->SetPoint( relIdx, (relIdx - 512)/SAMPLING_HELZ*1e6, a_hg_adc.at(j).at(hg_a_mainpeak_idx - 512 + relIdx )-pedestal_a_hg[j]  - (31-j)*100 );
            //         pLgGraph->SetPoint( relIdx, (relIdx - 512)/SAMPLING_HELZ*1e6, a_lg_adc.at(j).at(lg_a_mainpeak_idx - 512 + relIdx )-pedestal_a_lg[j]  - (31-j)*10  );
                    
            //         histHG.Fill( (relIdx - 512)/SAMPLING_HELZ*1e6 + 1, a_hg_adc.at(j).at(hg_a_mainpeak_idx - 512 + relIdx )-pedestal_a_hg[j] );
            //         histLG.Fill( (relIdx - 512)/SAMPLING_HELZ*1e6 + 1, a_lg_adc.at(j).at(lg_a_mainpeak_idx - 512 + relIdx )-pedestal_a_lg[j] );

            //     }
            //     hg_minoAnaWfTable.insert( std::make_pair( j, pHgGraph ) );
            //     lg_minoAnaWfTable.insert( std::make_pair( j, pLgGraph ) );
            //     h_all_mino_wf_hg->Add( &histHG );
            //     h_all_mino_wf_lg->Add( &histLG );
            // }
            std::cout << "test111" << std::endl;

			s->Delete();
			h_ROI->Delete();
		}
        std::cout << "test" << std::endl;

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
	g_a_hg_main_peak->SetMarkerSize(0.6);
	g_a_hg_main_peak->SetMarkerColor(kMagenta);
	g_a_main_rise->SetMarkerStyle(8);
	g_a_main_rise->SetMarkerSize(0.3);
	g_a_main_rise->SetMarkerColor(kRed);
	g_a_main_fall->SetMarkerStyle(8);
	g_a_main_fall->SetMarkerSize(0.3);
	g_a_main_fall->SetMarkerColor(kGreen);
	g_mino_peak->SetMarkerStyle(8);
	g_mino_peak->SetMarkerSize(0.3);
	g_mino_peak->SetMarkerColor(kRed);
	h_strip_a_hg->GetXaxis()->SetTitle("clock(1clock=0.4us)");
	h_strip_a_hg->GetYaxis()->SetTitle("CH");
	h_strip_a_hg->Draw("colz");
	g_a_hg_main_peak->Draw("same p");
	// g_a_main_rise->Draw("same p");
	// g_a_main_fall->Draw("same p");
	g_mino_peak->Draw("same p");
	
	c_strip->cd(2);
	g_a_lg_main_peak->SetMarkerStyle(8);
	g_a_lg_main_peak->SetMarkerSize(0.3);
	g_a_lg_main_peak->SetMarkerColor(kMagenta);
	h_strip_a_lg->GetXaxis()->SetTitle("clock(1clock=0.4us)");
	h_strip_a_lg->GetYaxis()->SetTitle("CH");
	h_strip_a_lg->GetZaxis()->SetRangeUser(-10.0, 10.0);
	h_strip_a_lg->Draw("colz");
	g_a_lg_main_peak->Draw("same p");
	
	c_strip->cd(3);
	g_c_hg_main_peak->SetMarkerStyle(8);
	g_c_hg_main_peak->SetMarkerSize(0.6);
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
	g_c_hg_main_peak->Draw("same p");
	g_c_main_rise->Draw("same p");
	g_c_main_fall->Draw("same p");
	
	c_strip->cd(4);
	g_c_lg_main_peak->SetMarkerStyle(8);
	g_c_lg_main_peak->SetMarkerSize(0.3);
	g_c_lg_main_peak->SetMarkerColor(kMagenta);
	h_strip_c_lg->GetXaxis()->SetTitle("clock(1clock=0.4us)");
	h_strip_c_lg->GetYaxis()->SetTitle("CH");
	h_strip_c_lg->Draw("colz");
	g_c_lg_main_peak->Draw("same p");

    c_strip->SaveAs( "recoil.png" );
    // c_strip->SaveAs( "recoil.pdf" );
    // c_strip->SaveAs( Form( "output6/strip_%d.png", ev_num ) );
	
	TCanvas *c_wave=new TCanvas("c_wave","",1000,1000);
	c_wave->Divide(2,2);
	// c_wave->cd(1)->DrawFrame(0,-4000,1600,4000,"HG waveform;us(4000sampling);mV");
	// c_wave->cd(2)->DrawFrame(0,-4000,1600,4000,"HG waveform;us(4000sampling);mV");

    TH1* pTmpHist = nullptr;
	// pTmpHist = c_wave->cd(1)->DrawFrame(900,-4000,1600,4000,"HG waveform;Time [us];Voltage (+offset) [mV]");
    // pTmpHist->GetYaxis()->SetTitleOffset( 1.7 );

	// pTmpHist = c_wave->cd(2)->DrawFrame(900,-4000,1600,4000,"HG waveform;Time [us];Voltage (+offset) [mV]");
    // pTmpHist->GetYaxis()->SetTitleOffset( 1.7 );

	pTmpHist = c_wave->cd(1)->DrawFrame(0,-3700,1600,3500,"HG waveform;Time [#mus];Amplitude [mV]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.7 );

	pTmpHist = c_wave->cd(2)->DrawFrame(0,-3700,1600,3500,"HG waveform;Time [#mus];Amplitude [mV]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.7 );

	pTmpHist = c_wave->cd(3)->DrawFrame(0,-350,1600,350,"LG waveform;Time [#mus];Amplitude [mV]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.5 );

	pTmpHist = c_wave->cd(4)->DrawFrame(0,-350,1600,350,"LG waveform;Time [#mus];Amplitude [mV]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.5 );

	for(int i=0;i<N_CHANNEL;i++){
		// c_wave->cd(1); if( i < 32 && mask_a_hg[i] == false ) g_a_hg[i]->Draw("p same");
		// c_wave->cd(3); if( i < 32 && mask_a_lg[i] == false ) g_a_lg[i]->Draw("p same");
		c_wave->cd(1); if( mask_a_hg[i] == false ) g_a_hg[i]->Draw("p same");
		c_wave->cd(3); if( mask_a_lg[i] == false ) g_a_lg[i]->Draw("p same");
		c_wave->cd(2); if( mask_c_hg[i] == false ) g_c_hg[i]->Draw("p same");
		c_wave->cd(4); if( mask_c_lg[i] == false ) g_c_lg[i]->Draw("p same");
	}

	c_wave->cd(1);
	g_mino_search->SetMarkerStyle(20);
	g_mino_search->SetMarkerSize(0.7);
	g_mino_search->SetMarkerColor(kGreen);
	g_mino_search->Draw("same p");

	g_mino_search_max->SetMarkerStyle(20);
	g_mino_search_max->SetMarkerSize(0.7);
	g_mino_search_max->SetMarkerColor(kRed);
	g_mino_search_max->Draw("same p");


    c_wave->SaveAs( "wave.png" );
    c_wave->SaveAs( "wave.pdf" );

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

    c_charge->SaveAs( "charge.png" );
    c_charge->SaveAs( "charge.pdf" );

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
	// app.Run();

    TCanvas* c_ave = new TCanvas( "c_ave", "c_ave", 800, 400 );
    c_ave->Divide( 2, 1 );
    
    pTmpHist = c_ave->cd(1)->DrawFrame(0, -47, 1600, 29, "Anode average waveform;Time [#mus];Amplitude [mV]" );
    pTmpHist->GetYaxis( )->SetTitleOffset( 1.5 );
    
    g_a_lg_ave->SetMarkerColor( kBlue );
    g_a_lg_ave->Draw( "p" );
    
    CreateDrawText( 0.17, 0.85, "Averaged waveform" );
    CreateDrawText( 0.17, 0.80, "Anode" );

    TPad* pSubPad1 = new TPad( "subpad1", "subpad1", 0.15, 0.15, 0.65, 0.6 );
    pSubPad1->Draw( );
    pTmpHist = pSubPad1->cd()->DrawFrame(825, -1, 1075, 27, "Anode average waveform;;" );
    pTmpHist->GetXaxis( )->SetLabelSize(0.07);
    pTmpHist->GetYaxis( )->SetLabelSize(0.07);
    TGraph* g_a_lg_ave_sub = dynamic_cast< TGraph* >( g_a_lg_ave->Clone( "a_lg_ave_sub" ) );
    g_a_lg_ave_sub->Draw( "p" );

    pTmpHist = c_ave->cd(2)->DrawFrame(0, -47, 1600, 29, "Cathode average waveform;Time [#mus];Amplitude [mV]" );
    pTmpHist->GetYaxis( )->SetTitleOffset( 1.5 );
    g_c_lg_ave->SetMarkerColor( kBlue );
    g_c_lg_ave->Draw( "p" );

    CreateDrawText( 0.17, 0.85, "Averaged waveform" );
    CreateDrawText( 0.17, 0.80, "Cathode" );

    TPad* pSubPad2 = new TPad( "subpad2", "subpad2", 0.15, 0.15, 0.65, 0.6 );
    pSubPad2->Draw( );
    pTmpHist = pSubPad2->cd()->DrawFrame(825, -27, 1075, 1, "Anode average waveform;;" );
    pTmpHist->GetXaxis( )->SetLabelSize(0.07);
    pTmpHist->GetYaxis( )->SetLabelSize(0.07);
    TGraph* g_c_lg_ave_sub = dynamic_cast< TGraph* >( g_c_lg_ave->Clone( "c_lg_ave_sub" ) );
    g_c_lg_ave_sub->Draw( "p" );

    c_ave->SaveAs( "ave_lg.png" );
    c_ave->SaveAs( "ave_lg.pdf" );


    pTmpHist = c_ave->cd(1)->DrawFrame(0, -1000, 1600, 600, "Anode average waveform;Time [#mus];Amplitude [mV]" );
    pTmpHist->GetYaxis( )->SetTitleOffset( 1.5 );
    
    g_a_hg_ave->SetMarkerColor( kBlue );
    g_a_hg_ave->Draw( "p" );
    
    CreateDrawText( 0.17, 0.85, "Averaged waveform" );
    CreateDrawText( 0.17, 0.80, "Anode" );

    TPad* pSubPad1hg = new TPad( "subpad1hg", "subpad1hg", 0.15, 0.15, 0.65, 0.6 );
    pSubPad1hg->Draw( );
    pTmpHist = pSubPad1hg->cd()->DrawFrame(1025, -20, 1275, 200, "Anode average waveform;;" );
    pTmpHist->GetXaxis( )->SetLabelSize(0.07);
    pTmpHist->GetYaxis( )->SetLabelSize(0.07);
    TGraph* g_a_hg_ave_sub = dynamic_cast< TGraph* >( g_a_hg_ave->Clone( "a_hg_ave_sub" ) );
    g_a_hg_ave_sub->Draw( "p" );

    pTmpHist = c_ave->cd(2)->DrawFrame(0, -1000, 1600, 600, "Cathode average waveform;Time [#mus];Amplitude [mV]" );
    pTmpHist->GetYaxis( )->SetTitleOffset( 1.5 );
    g_c_hg_ave->SetMarkerColor( kBlue );
    g_c_hg_ave->Draw( "p" );

    CreateDrawText( 0.17, 0.85, "Averaged waveform" );
    CreateDrawText( 0.17, 0.80, "Cathode" );

    TPad* pSubPad2hg = new TPad( "subpad2hg", "subpad2hg", 0.15, 0.15, 0.65, 0.6 );
    pSubPad2hg->Draw( );
    pTmpHist = pSubPad2hg->cd()->DrawFrame(825, -600, 1075, 20, "Anode average waveform;;" );
    pTmpHist->GetXaxis( )->SetLabelSize(0.07);
    pTmpHist->GetYaxis( )->SetLabelSize(0.07);
    TGraph* g_c_hg_ave_sub = dynamic_cast< TGraph* >( g_c_hg_ave->Clone( "c_hg_ave_sub" ) );
    g_c_hg_ave_sub->Draw( "p" );

    c_ave->SaveAs( "ave_hg.png" );
    c_ave->SaveAs( "ave_hg.pdf" );



    TCanvas* c_image = new TCanvas( "c_image", "c_image", 800, 400 );
    c_image->Divide( 2, 1 );
    c_image->cd( 1 );
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.80, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
    gStyle->SetNumberContours( NCont );
    gPad->SetRightMargin( 0.2 );
	h_track_a_lg->SetXTitle( "Time [#mus]" );
	h_track_a_lg->SetYTitle( "x [cm]" );
	h_track_a_lg->SetZTitle( "Amplitude [mV]" );
	h_track_a_lg->GetYaxis()->SetRangeUser(0, 1.28);
	h_track_a_lg->GetZaxis()->SetTitleOffset( 1.4 );
	h_track_a_lg->Draw("colz");

    c_image->cd( 2 );
    TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
    gStyle->SetNumberContours( NCont );
    gPad->SetRightMargin( 0.2 );
	h_track_c_lg->SetXTitle( "Time [#mus]" );
	h_track_c_lg->SetYTitle( "y [cm]" );
	h_track_c_lg->SetZTitle( "Amplitude [mV]" );
	h_track_c_lg->GetZaxis()->SetTitleOffset( 1.4 );
	h_track_c_lg->Draw("colz");

    c_image->SaveAs( "image.png" );
    c_image->SaveAs( "image.pdf" );

    TCanvas* c_dt = new TCanvas( "c_dt", "c_dt", 800, 600 );
    h_dt->SetXTitle( "#Deltat [#mus]" );
    h_dt->SetYTitle( "Entries" );
    h_dt->Draw( );
    TF1 fitGauss( "fitGauss", "gaus", 0, 200 );
    fitGauss.SetParameter( 1, h_dt->GetMean( ) );
    fitGauss.SetParameter( 2, h_dt->GetRMS( ) );
    fitGauss.SetParLimits( 1, 0.0, 200.0 );
    h_dt->Fit( &fitGauss, "LM", "" );

    CreateDrawText( 0.6, 0.88, Form( "#Deltat = %0.2lf", fitGauss.GetParameter( 1 ) ) );
    CreateDrawText( 0.6, 0.8, Form( "#sigma_{#Deltat} = %0.2lf", fitGauss.GetParameter( 2 ) ) );
    CreateDrawText( 0.6, 0.72, Form( "chi2 / NDF = %0.2lf / %d", fitGauss.GetChisquare( ), fitGauss.GetNDF( ) ) );

    c_dt->SaveAs( "dt_1.png" );
    c_dt->SaveAs( "dt_1.pdf" );

    h_dt->Draw( );
    fitGauss.SetParLimits( 1, fitGauss.GetParameter( 1 ) - fitGauss.GetParameter( 2 )*1.2, fitGauss.GetParameter( 1 ) + fitGauss.GetParameter( 2 )*1.2 );
    h_dt->Fit( &fitGauss, "LM", "", fitGauss.GetParameter( 1 ) - fitGauss.GetParameter( 2 )*1.2, fitGauss.GetParameter( 1 ) + fitGauss.GetParameter( 2 )*1.2 );

    CreateDrawText( 0.6, 0.88, Form( "#Deltat = %0.2lf", fitGauss.GetParameter( 1 ) ) );
    CreateDrawText( 0.6, 0.8, Form( "#sigma_{#Deltat} = %0.2lf", fitGauss.GetParameter( 2 ) ) );
    CreateDrawText( 0.6, 0.72, Form( "chi2 / NDF = %0.2lf / %d", fitGauss.GetChisquare( ), fitGauss.GetNDF( ) ) );

    c_dt->SaveAs( "dt_2.png" );
    c_dt->SaveAs( "dt_2.pdf" );

    h_dt->Draw( );
    fitGauss.SetParLimits( 1, fitGauss.GetParameter( 1 ) - fitGauss.GetParameter( 2 ), fitGauss.GetParameter( 1 ) + fitGauss.GetParameter( 2 ) );
    h_dt->Fit( &fitGauss, "LM", "", fitGauss.GetParameter( 1 ) - fitGauss.GetParameter( 2 ), fitGauss.GetParameter( 1 ) + fitGauss.GetParameter( 2 ) );

    CreateDrawText( 0.6, 0.88, Form( "#Deltat = %0.2lf", fitGauss.GetParameter( 1 ) ) );
    CreateDrawText( 0.6, 0.8, Form( "#sigma_{#Deltat} = %0.2lf", fitGauss.GetParameter( 2 ) ) );
    CreateDrawText( 0.6, 0.72, Form( "chi2 / NDF = %0.2lf / %d", fitGauss.GetChisquare( ), fitGauss.GetNDF( ) ) );

    c_dt->SaveAs( "dt_3.png" );
    c_dt->SaveAs( "dt_3.pdf" );


    TCanvas* c_minoAnaWf = new TCanvas( "c_minoAnaWf", "c_minoAnaWf", 800, 400 );
    c_minoAnaWf->Divide( 2, 1 );

	pTmpHist = c_minoAnaWf->cd(1)->DrawFrame(-200,-3000,200,1000,"HG waveform;Relative time from main peak [#mus];Amplitude [mV]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.7 );
    for( auto pair : hg_minoAnaWfTable ) {
        if( pair.second != nullptr ) {
            pair.second->SetMarkerColor( MyPalette[pair.first] );
            pair.second->SetMarkerStyle( 1 );
            pair.second->Draw( "same p" );
        }
    }

	pTmpHist = c_minoAnaWf->cd(2)->DrawFrame(-200,-300,200,50,"LG waveform;Relative time from main peak [#mus];Amplitude [mV]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.5 );
    for( auto pair : lg_minoAnaWfTable ) {
        if( pair.second != nullptr ) {
            pair.second->SetMarkerColor( MyPalette[pair.first] );
            pair.second->SetMarkerStyle( 1 );
            pair.second->Draw( "same p" );
        }
    }

    c_minoAnaWf->SaveAs( "minoAna.png" );
    c_minoAnaWf->SaveAs( "minoAna.pdf" );


    h_all_mino_wf_hg->Scale( 1.0 / h_all_mino_wf_hg->GetMaximum( ) );
    h_all_mino_wf_lg->Scale( 1.0 / h_all_mino_wf_lg->GetMaximum( ) );

	pTmpHist = c_minoAnaWf->cd(1)->DrawFrame(-200,-0.2,200,1.2,"HG waveform;Relative time from main peak [#mus];Summed amplitude [A.U.]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.7 );
    h_all_mino_wf_hg->SetMarkerStyle( 8 );
    h_all_mino_wf_hg->SetMarkerSize( 0.5 );
    h_all_mino_wf_hg->SetMarkerColor( kBlue );
    h_all_mino_wf_hg->Draw( "same histp" );

	pTmpHist = c_minoAnaWf->cd(2)->DrawFrame(-200,-0.2,200,1.2,"LG waveform;Relative time from main peak [#mus];Summed amplitude [A.U.]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.5 );
    h_all_mino_wf_lg->SetMarkerStyle( 8 );
    h_all_mino_wf_lg->SetMarkerSize( 0.5 );
    h_all_mino_wf_lg->SetMarkerColor( kBlue );
    h_all_mino_wf_lg->Draw( "same histp" );
    
    c_minoAnaWf->SaveAs( "allMinoAna.png" );
    c_minoAnaWf->SaveAs( "allMinoAna.pdf" );


	pTmpHist = c_minoAnaWf->cd(1)->DrawFrame(-200,-0.1,50,0.3,"HG waveform;Relative time from main peak [#mus];Summed amplitude [A.U.]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.7 );
    h_all_mino_wf_hg->SetMarkerStyle( 8 );
    h_all_mino_wf_hg->SetMarkerSize( 0.5 );
    h_all_mino_wf_hg->SetMarkerColor( kBlue );
    h_all_mino_wf_hg->Draw( "same histp" );

	pTmpHist = c_minoAnaWf->cd(2)->DrawFrame(-200,-0.1,50,0.3,"LG waveform;Relative time from main peak [#mus];Summed amplitude [A.U.]");
    pTmpHist->GetYaxis()->SetTitleOffset( 1.5 );
    h_all_mino_wf_lg->SetMarkerStyle( 8 );
    h_all_mino_wf_lg->SetMarkerSize( 0.5 );
    h_all_mino_wf_lg->SetMarkerColor( kBlue );
    h_all_mino_wf_lg->Draw( "same histp" );
    
    c_minoAnaWf->SaveAs( "allMinoAnaZoom.png" );
    c_minoAnaWf->SaveAs( "allMinoAnaZoom.pdf" );


	return 1;
}
