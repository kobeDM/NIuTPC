
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
#include <numeric>
// ROOT
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
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
	Double_t Red[]    = {1.0, 0.0, 1.0};
	Double_t Green[]  = {0.,  0.0, 0.};
	Double_t Blue[]   = {0.,  1.0, 1.0};
	Double_t Length[] = {0., .50, 1.0};
	Int_t FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 64);
	for (int i=0;i<64;i++) MyPalette[i] = FI+i;

	if(argc <3){
		std::cerr << "Usage:" << std::endl;
		std::cerr << "/.main [filename.root] [config.json]" << std::endl;
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
	int offset_sampling       = ni_conf->offset_sampling;
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

	//++++++++++++++++++++++++++++++++++++++++++
	//  Iutput Tree
	//++++++++++++++++++++++++++++++++++++++++++
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


	//++++++++++++++++++++++++++++++++++++++++++
	//  Output Tree
	//++++++++++++++++++++++++++++++++++++++++++
	std::string filename;
	std::string::size_type pos = dirfilename.find("/");
	while(pos != std::string::npos){
		dirfilename = dirfilename.substr(pos+1);
		pos = dirfilename.find("/");
	}
	pos = dirfilename.find(".root");
	filename = dirfilename.substr(0,pos);
	std::string outfilename = filename.substr(0, index) + "_anal.root";
	TFile *fout=new TFile(outfilename.c_str(),"recreate");
	TTree* outtree = new TTree("anal_tree","anal_tree");
	int outtree_ev_num;
	//anode
	double outtree_a_hg_sum_pulse_height,outtree_a_lg_sum_pulse_height;
	double outtree_a_hg_sum_charge,outtree_a_lg_sum_charge;
	double outtree_a_hg_sum_tot;
	vector<double> outtree_a_hg_pulse_height,outtree_a_lg_pulse_height;
	vector<double> outtree_a_hg_charge,outtree_a_lg_charge;
	vector<double> outtree_a_hg_mainpeak_time,outtree_a_lg_mainpeak_time;
	vector<double> outtree_a_hg_mainrise_time;
	vector<double> outtree_a_hg_mainfall_time;
	vector<double> outtree_a_hg_tot;
	//cathode
	double outtree_c_hg_sum_pulse_height,outtree_c_lg_sum_pulse_height;
	double outtree_c_hg_sum_charge,outtree_c_lg_sum_charge;
	double outtree_c_hg_sum_tot;
	vector<double> outtree_c_hg_pulse_height,outtree_c_lg_pulse_height;
	vector<double> outtree_c_hg_charge,outtree_c_lg_charge;
	vector<double> outtree_c_hg_mainpeak_time,outtree_c_lg_mainpeak_time;
	vector<double> outtree_c_hg_mainrise_time;
	vector<double> outtree_c_hg_mainfall_time;
	vector<double> outtree_c_hg_tot;
	//detector
	vector<double> outtree_xz_x,outtree_xz_z,outtree_yz_y,outtree_yz_z;
	double outtree_ave_x,outtree_ave_y,outtree_ave_z;
	//minority (absolute z reconstruction)
	vector<double> outtree_a_hg_minopeak_time;
	vector<double> outtree_c_hg_minopeak_time;
	vector<double> outtree_dt;
	vector<double> outtree_abs_x,outtree_abs_y,outtree_abs_z;
	double outtree_ave_abs_x, outtree_ave_abs_y, outtree_ave_abs_z;

	outtree->Branch("ev",                   &outtree_ev_num);
	outtree->Branch("a_hg_sum_pulse_height",&outtree_a_hg_sum_pulse_height);
	outtree->Branch("a_lg_sum_pulse_height",&outtree_a_lg_sum_pulse_height);
	outtree->Branch("a_hg_sum_charge",      &outtree_a_hg_sum_charge);
	outtree->Branch("a_lg_sum_charge",      &outtree_a_lg_sum_charge);
	outtree->Branch("a_hg_sum_tot",         &outtree_a_hg_sum_tot);
	outtree->Branch("a_hg_pulse_height",    &outtree_a_hg_pulse_height);
	outtree->Branch("a_lg_pulse_height",    &outtree_a_lg_pulse_height);
	outtree->Branch("a_hg_charge",          &outtree_a_hg_charge);
	outtree->Branch("a_lg_charge",          &outtree_a_lg_charge);
	outtree->Branch("a_hg_mainpeak_time",   &outtree_a_hg_mainpeak_time);
	outtree->Branch("a_lg_mainpeak_time",   &outtree_a_lg_mainpeak_time);
	outtree->Branch("a_hg_mainrise_time",   &outtree_a_hg_mainrise_time);
	outtree->Branch("a_hg_mainfall_time",   &outtree_a_hg_mainfall_time);
	outtree->Branch("a_hg_tot",             &outtree_a_hg_tot);

	outtree->Branch("c_hg_sum_pulse_height",&outtree_c_hg_sum_pulse_height);
	outtree->Branch("c_lg_sum_pulse_height",&outtree_c_lg_sum_pulse_height);
	outtree->Branch("c_hg_sum_charge",      &outtree_c_hg_sum_charge);
	outtree->Branch("c_lg_sum_charge",      &outtree_c_lg_sum_charge);
	outtree->Branch("c_hg_sum_tot",         &outtree_c_hg_sum_tot);
	outtree->Branch("c_hg_pulse_height",    &outtree_c_hg_pulse_height);
	outtree->Branch("c_lg_pulse_height",    &outtree_c_lg_pulse_height);
	outtree->Branch("c_hg_charge",          &outtree_c_hg_charge);
	outtree->Branch("c_lg_charge",          &outtree_c_lg_charge);
	outtree->Branch("c_hg_mainpeak_time",   &outtree_c_hg_mainpeak_time);
	outtree->Branch("c_lg_mainpeak_time",   &outtree_c_lg_mainpeak_time);
	outtree->Branch("c_hg_mainrise_time",   &outtree_c_hg_mainrise_time);
	outtree->Branch("c_hg_mainfall_time",   &outtree_c_hg_mainfall_time);
	outtree->Branch("c_hg_tot",             &outtree_c_hg_tot);

	outtree->Branch("xz_x",                 &outtree_xz_x);
	outtree->Branch("xz_z",                 &outtree_xz_z);
	outtree->Branch("yz_y",                 &outtree_yz_y);
	outtree->Branch("yz_z",                 &outtree_yz_z);
	outtree->Branch("ave_x",                &outtree_ave_x);
	outtree->Branch("ave_y",                &outtree_ave_y);
	outtree->Branch("ave_z",                &outtree_ave_z);

	outtree->Branch("a_hg_minopeak_time",   &outtree_a_hg_minopeak_time);
	outtree->Branch("c_hg_minopeak_time",   &outtree_c_hg_minopeak_time);
	outtree->Branch("dt",                   &outtree_dt);
	outtree->Branch("abs_x",                &outtree_abs_x);
	outtree->Branch("abs_y",                &outtree_abs_y);
	outtree->Branch("abs_z",                &outtree_abs_z);
	outtree->Branch("ave_abs_x",            &outtree_ave_abs_x);
	outtree->Branch("ave_abs_y",            &outtree_ave_abs_y);
	outtree->Branch("ave_abs_z",            &outtree_ave_abs_z);

	
	//++++++++++++++++++++++++++++++++++++++++++
	//  visualize
	//++++++++++++++++++++++++++++++++++++++++++
	TH1D* h_sum_pulse_height = new TH1D("h_sum_pulse_height","h_sum_pulse_height",40,0,5000);
	TH1D* h_sum_charge       = new TH1D("h_sum_charge"      ,"h_sum_charge"      ,40,0,250000);
	TH2D* h_xy               = new TH2D("h_xy"              ,"h_xy"              ,64,-1.4,1.16,64,-1.4,1.16); //bin width = 400um
	TH2D* h_xz               = new TH2D("h_xz"              ,"h_xz"              ,64,-1.4,1.16,360,0,14.4); //bin width = 400um
	TH2D* h_yz               = new TH2D("h_yz"              ,"h_yz"              ,360,0,14.4,64,-1.4,1.16); //bin width = 400um
	TH1D* h_dt               = new TH1D("h_dt"              ,"h_dt"              ,1000,0,400); //bin width 0.4us
	TH1D* h_abs_z            = new TH1D("h_abs_z"           ,"h_abs_z"           ,64,0,14.4); //bin width = 400um
	TH1D* h_strip_dt         = new TH1D("h_strip_dt"        ,"h_strip_dt"        ,1000,0,400); //bin width 0.4us
	TH1D* h_strip_abs_z      = new TH1D("h_strip_abs_z"     ,"h_strip_abs_z"     ,64,0,14.4); //bin width = 400um

	//for debug
	TH1D* h_31ch_32ch_dt     = new TH1D("h_31ch_32ch_dt"    ,"h_31ch_32ch_dt"    ,SAMPLING_NUM,-SAMPLING_NUM,SAMPLING_NUM);

	//++++++++++++++++++++++++++++++++++++++++++
	//  main loop
	//++++++++++++++++++++++++++++++++++++++++++
	for(int ev=0;ev<nevent;ev++){
	//for(int ev=0;ev<1;ev++){ //to test
		// Get Event Info
		tree->GetEntry(ev);
		if(ev%1==0) std::cerr << "\rFill Data into Tree ... : "<< ev << "/" << nevent << std::flush;
		
		//vector initialization
		outtree_a_hg_pulse_height.clear();
		outtree_a_lg_pulse_height.clear();
		outtree_a_hg_charge.clear();
		outtree_a_lg_charge.clear();
		outtree_a_hg_mainpeak_time.clear();
		outtree_a_lg_mainpeak_time.clear();
		outtree_a_hg_mainrise_time.clear();
		outtree_a_hg_mainfall_time.clear();
		outtree_c_hg_pulse_height.clear();
		outtree_c_lg_pulse_height.clear();
		outtree_c_hg_charge.clear();
		outtree_c_lg_charge.clear();
		outtree_c_hg_mainpeak_time.clear();
		outtree_c_lg_mainpeak_time.clear();
		outtree_c_hg_mainrise_time.clear();
		outtree_c_hg_mainfall_time.clear();
		outtree_xz_x.clear();
		outtree_xz_z.clear();
		outtree_yz_y.clear();
		outtree_yz_z.clear();
        outtree_a_hg_minopeak_time.clear();
		outtree_c_hg_minopeak_time.clear();
		outtree_dt.clear();
		outtree_abs_x.clear();
		outtree_abs_y.clear();
		outtree_abs_z.clear();

		int trigger_num = event->GetTrigger();

        outtree_ev_num = ev;
		vector<vector<double>> a_hg_adc = event->GetAnodeHGADC();
		vector<vector<double>> a_lg_adc = event->GetAnodeLGADC();
		vector<vector<double>> c_hg_adc = event->GetCathodeHGADC();
		vector<vector<double>> c_lg_adc = event->GetCathodeLGADC();

		// Pedestal calc
		double pedestal_a_hg[64];
		double pedestal_a_lg[64];
		double pedestal_c_hg[64];
		double pedestal_c_lg[64];
		for(int j=0;j<N_CHANNEL;j++){
			pedestal_a_hg[j]=0;
			pedestal_a_lg[j]=0;
			pedestal_c_hg[j]=0;
			pedestal_c_lg[j]=0;
			for(int k=0;k<offset_sampling;k++){
				pedestal_a_hg[j] += a_hg_adc[j][k]/offset_sampling;
				pedestal_a_lg[j] += a_lg_adc[j][k]/offset_sampling;
				pedestal_c_hg[j] += c_hg_adc[j][k]/offset_sampling;
				pedestal_c_lg[j] += c_lg_adc[j][k]/offset_sampling;
			}
		}

		//anode
		double sum_a_hg_pulse_height = 0;
		double sum_a_lg_pulse_height = 0;
		double sum_a_hg_charge = 0;
		double sum_a_lg_charge = 0;
		double sum_a_hg_tot = 0;
		//cathode
		double sum_c_hg_pulse_height = 0;
		double sum_c_lg_pulse_height = 0;
		double sum_c_hg_charge = 0;
		double sum_c_lg_charge = 0;
		double sum_c_hg_tot = 0;

		//for debug
		double peak_time_31ch = sqrt(-1);;
		double peak_time_32ch = sqrt(-1);;

		for(int j=0;j<N_CHANNEL;j++){
			//++++++++++++++++++++++++++++++++++++++++++
			//  Main Peak analysis
			//++++++++++++++++++++++++++++++++++++++++++
			//anode
			double hg_a_mainpeak_time=-1;
			double lg_a_mainpeak_time=-1;
			double hg_a_mainrise_time=-1;
			double lg_a_mainrise_time=-1;
			double hg_a_mainfall_time=-1;
			double lg_a_mainfall_time=-1;
			double hg_a_pulse_max=0;
			double lg_a_pulse_max=0;
			double hg_a_this_charge = 0;
			double lg_a_this_charge = 0;
			double hg_a_this_tot = 0;
			int hg_a_rise_flag=0;
			//cathode
			double hg_c_mainpeak_time=-1;
			double lg_c_mainpeak_time=-1;
			double hg_c_mainrise_time=-1;
			double lg_c_mainrise_time=-1;
			double hg_c_mainfall_time=-1;
			double lg_c_mainfall_time=-1;
			double hg_c_pulse_max=0;
			double lg_c_pulse_max=0;
			double hg_c_this_charge = 0;
			double lg_c_this_charge = 0;
			double hg_c_this_tot = 0;
			int hg_c_rise_flag=0;
			//for minority
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
					hg_a_this_tot += 1./SAMPLING_HELZ*1e6; //us
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
					hg_c_this_tot += 1./SAMPLING_HELZ*1e6; //us
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
					//maximum method
					if(some_mino_height[pks]>mino_peak_max){
						this_mino_time = some_mino_time[pks];
						mino_peak_max = some_mino_height[pks];
					}
					/*
					// dt minimum method
					if(lg_a_mainpeak_time-some_mino_time[pks]<mino_dt_min){
						this_mino_time = some_mino_time[pks];
						mino_dt_min = lg_a_mainpeak_time-some_mino_time[pks]; 
					}
					*/
				}

				s->Delete();
				h_ROI->Delete();
			}
			h_mino_search->Delete();
			
			//++++++++++++++++++++++++++++++++++++++++++
			//  Position analysis
			//  data push_back
			//++++++++++++++++++++++++++++++++++++++++++
			//main peak triggered
			if(hg_a_mainpeak_time!=-1){
				double this_xz_x = (29-j+0.5)*0.04;//cm
				double this_xz_z = lg_a_mainpeak_time*driftV_main*1e-3;//cm
				//std::cout << "x : " << this_xz_x << " cm , z:" << this_xz_z << " cm" <<std::endl;
				//anode
				outtree_a_hg_pulse_height.push_back(hg_a_pulse_max);
				outtree_a_lg_pulse_height.push_back(lg_a_pulse_max);
				outtree_a_hg_charge.push_back(hg_a_this_charge);
				outtree_a_lg_charge.push_back(lg_a_this_charge);
				sum_a_hg_pulse_height +=hg_a_pulse_max;
				sum_a_lg_pulse_height +=lg_a_pulse_max;
				sum_a_hg_charge +=hg_a_this_charge;
				sum_a_lg_charge +=lg_a_this_charge;
				sum_a_hg_tot    +=hg_a_this_tot;
				outtree_a_hg_mainpeak_time.push_back(hg_a_mainpeak_time);
				outtree_a_lg_mainpeak_time.push_back(lg_a_mainpeak_time);
				outtree_a_hg_mainrise_time.push_back(hg_a_mainrise_time);
				outtree_a_hg_mainfall_time.push_back(hg_a_mainfall_time);
				outtree_xz_x.push_back(this_xz_x);
				outtree_xz_z.push_back(this_xz_z);
				outtree_a_hg_tot.push_back(hg_a_this_tot);
				h_xz->Fill(this_xz_x,this_xz_z);
				//minority triggered
				if(this_mino_time!=-1){
					//if(hg_a_pulse_max<300)continue;
					double this_abs_z = (lg_a_mainpeak_time-this_mino_time)*calc_abs_z_param*1e-3;
					h_strip_dt->Fill(lg_a_mainpeak_time-this_mino_time);
					h_strip_abs_z->Fill(this_abs_z);
					//std::cout << "delta t : " << lg_a_mainpeak_time-this_mino_time << " us , calc param:" << calc_abs_z_param << "" <<std::endl;
					//std::cout << "abs x : " << this_xz_x << " cm , z:" << this_abs_z << " cm" <<std::endl;
					outtree_dt.push_back(lg_a_mainpeak_time-this_mino_time);
					outtree_abs_x.push_back(this_xz_x);
					outtree_abs_z.push_back(this_abs_z);
					outtree_a_hg_minopeak_time.push_back(this_mino_time);
				}
			}
			if(hg_c_mainpeak_time!=-1){
				double this_yz_y = (29-j+0.5)*0.04;//cm
				double this_yz_z = lg_c_mainpeak_time*driftV_main*1e-3;//cm
				//std::cout << "y : " << this_yz_y << " cm , z:" << this_yz_z << " cm" <<std::endl;
				//cathode
				outtree_c_hg_pulse_height.push_back(hg_c_pulse_max);
				outtree_c_lg_pulse_height.push_back(lg_c_pulse_max);
				outtree_c_hg_charge.push_back(hg_c_this_charge);
				outtree_c_lg_charge.push_back(lg_c_this_charge);
				sum_c_hg_pulse_height +=hg_c_pulse_max;
				sum_c_lg_pulse_height +=lg_c_pulse_max;
				sum_c_hg_charge +=hg_c_this_charge;
				sum_c_lg_charge +=lg_c_this_charge;
				sum_c_hg_tot    +=hg_c_this_tot;
				outtree_c_hg_mainpeak_time.push_back(hg_c_mainpeak_time);
				outtree_c_lg_mainpeak_time.push_back(lg_c_mainpeak_time);
				outtree_c_hg_mainrise_time.push_back(hg_c_mainrise_time);
				outtree_c_hg_mainfall_time.push_back(hg_c_mainfall_time);
				outtree_yz_y.push_back(this_yz_y);
				outtree_yz_z.push_back(this_yz_z);
				outtree_c_hg_tot.push_back(hg_c_this_tot);
				h_yz->Fill(this_yz_z,this_yz_y);
			}
			if(j==31)peak_time_31ch = lg_a_mainpeak_time;
			if(j==32)peak_time_32ch = lg_a_mainpeak_time;
		}
		//++++++++++++++++++++++++++++++++++++++++++
		//  Fill Data
		//++++++++++++++++++++++++++++++++++++++++++
		
		outtree_a_hg_sum_pulse_height = sum_a_hg_pulse_height;
		outtree_a_lg_sum_pulse_height = sum_a_lg_pulse_height;
		outtree_a_hg_sum_charge = sum_a_hg_charge;
		outtree_a_lg_sum_charge = sum_a_lg_charge;
		outtree_a_hg_sum_tot = sum_a_hg_tot;
		outtree_c_hg_sum_pulse_height = sum_c_hg_pulse_height;
		outtree_c_lg_sum_pulse_height = sum_c_lg_pulse_height;
		outtree_c_hg_sum_charge = sum_c_hg_charge;
		outtree_c_lg_sum_charge = sum_c_lg_charge;
		outtree_c_hg_sum_tot = sum_c_hg_tot;

		if(outtree_xz_x.size()==0 && outtree_yz_y.size()==0){
			outtree_ave_x = double(sqrt(-1));
			outtree_ave_y = double(sqrt(-1));
			outtree_ave_z = double(sqrt(-1));
			continue;
		}else if(outtree_xz_x.size()==0){
			outtree_ave_x = double(sqrt(-1));
			outtree_ave_y = std::accumulate(outtree_yz_y.begin(),outtree_yz_y.end(),0.0)/outtree_yz_y.size();
			outtree_ave_z = std::accumulate(outtree_yz_z.begin(),outtree_yz_z.end(),0.0)/outtree_yz_z.size();
			//continue;
		}else if(outtree_yz_y.size()==0){
			outtree_ave_x = std::accumulate(outtree_xz_x.begin(),outtree_xz_x.end(),0.0)/outtree_xz_x.size();
			outtree_ave_y = double(sqrt(-1));
			outtree_ave_z = std::accumulate(outtree_xz_z.begin(),outtree_xz_z.end(),0.0)/outtree_xz_z.size();
			//continue;
		}else{
			outtree_ave_x = std::accumulate(outtree_xz_x.begin(),outtree_xz_x.end(),0.0)/outtree_xz_x.size();
			outtree_ave_y = std::accumulate(outtree_yz_y.begin(),outtree_yz_y.end(),0.0)/outtree_yz_y.size();
			outtree_ave_z = (std::accumulate(outtree_xz_z.begin(),outtree_xz_z.end(),0.0)/outtree_xz_z.size()
											+std::accumulate(outtree_yz_z.begin(),outtree_yz_z.end(),0.0)/outtree_yz_z.size()
											)/2;
			h_xy->Fill(outtree_ave_x,outtree_ave_y);
			//std::cout << "xz.size() : " << outtree_xz_x.size() << " | ave x : " << outtree_ave_x << "cm | ave z : " << outtree_ave_z << "cm" <<std::endl; 
			//std::cout << "ave y : " << outtree_ave_y << "cm | ave z : " << outtree_ave_z << "cm" <<std::endl; 
		}
		if(outtree_abs_z.size()==0){
			outtree_ave_abs_x = double(sqrt(-1));
			outtree_ave_abs_z = double(sqrt(-1));
		}else{
			outtree_ave_abs_x = std::accumulate(outtree_abs_x.begin(),outtree_abs_x.end(),0.0)/outtree_abs_x.size();
			outtree_ave_abs_z = std::accumulate(outtree_abs_z.begin(),outtree_abs_z.end(),0.0)/outtree_abs_z.size();
			//std::cout << "abs x : " << outtree_ave_abs_x << "cm | abs z : " << outtree_ave_abs_z << "cm" <<std::endl; 
			h_abs_z->Fill(outtree_ave_abs_z);
			h_dt->Fill(outtree_ave_abs_z*1e3/calc_abs_z_param);
		}
		outtree->Fill();
		h_sum_pulse_height->Fill(sum_a_lg_pulse_height);
		h_sum_charge->Fill(sum_a_lg_charge);
	
		//for debug
		h_31ch_32ch_dt->Fill(peak_time_32ch-peak_time_31ch);
        std::cout << dirfilename << "\t" << ev << std::endl;
	}

	//++++++++++++++++++++++++++++++++++++++++++
	//  Draw
	//++++++++++++++++++++++++++++++++++++++++++
	gStyle->SetPalette(kRainBow);
	TCanvas* c_vis = new TCanvas("c_vis","c_vis",0,0,1000,1000);
	c_vis->Divide(3,3);
	c_vis->cd(2);
	h_sum_pulse_height->GetXaxis()->SetTitle("sum_pulse_height(mV)");
	h_sum_pulse_height->GetYaxis()->SetTitle("counts");
	h_sum_pulse_height->Draw();
	c_vis->cd(5);
	h_sum_charge->GetXaxis()->SetTitle("sum_charge");
	h_sum_charge->GetYaxis()->SetTitle("counts");
	h_sum_charge->Draw();
	c_vis->cd(4);
	h_xz->GetXaxis()->SetTitle("x(cm)");
	h_xz->GetYaxis()->SetTitle("z(cm)");
	h_xz->Draw("colz");
	c_vis->cd(7);
	h_xy->GetXaxis()->SetTitle("x(cm)");
	h_xy->GetYaxis()->SetTitle("y(cm)");
	h_xy->Draw("colz");
	c_vis->cd(8);
	h_yz->GetXaxis()->SetTitle("z(cm)");
	h_yz->GetYaxis()->SetTitle("y(cm)");
	h_yz->Draw("colz");

	TCanvas* c_mino = new TCanvas("c_mino","c_mino",0,0,1000,1000);
	c_mino->Divide(3,3);
	c_mino->cd(5);
	h_strip_dt->GetXaxis()->SetTitle("main - minority(us)");
	h_strip_dt->GetYaxis()->SetTitle("counts");
	h_strip_dt->Draw();
	c_mino->cd(8);
	h_strip_abs_z->GetXaxis()->SetTitle("absolute z(cm)");
	h_strip_abs_z->GetYaxis()->SetTitle("counts");
	h_strip_abs_z->Draw();
	c_mino->cd(6);
	h_dt->GetXaxis()->SetTitle("main - minority(us)");
	h_dt->GetYaxis()->SetTitle("counts");
	h_dt->Draw();
	c_mino->cd(9);
	h_abs_z->GetXaxis()->SetTitle("absolute z(cm)");
	h_abs_z->GetYaxis()->SetTitle("counts");
	h_abs_z->Draw();

	TCanvas* c_debug = new TCanvas("c_debug","c_debug",0,0,1000,1000);
	h_31ch_32ch_dt->GetXaxis()->SetTitle("#Delta clock(2.5Mhz)");
	h_31ch_32ch_dt->GetYaxis()->SetTitle("counts");
	h_31ch_32ch_dt->Draw();

	outtree->Write();
	// app.Run();
	//++++++++++++++++++++++++++++++++++++++++++
	//  Output
	//++++++++++++++++++++++++++++++++++++++++++
	fout->Close();
	std::cout << "" <<std::endl;
	std::cout << "--- anal end ---" <<std::endl;

	// app.Run();

	return 1;
}




