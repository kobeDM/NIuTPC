
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
	std::string::size_type index = dirfilename.find("_anal.root");
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
	int offset_sampling          = ni_conf->offset_sampling;
	double cal_factor            = ni_conf->cal_factor;
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
	TApplication app("app",&argc,argv);  


	//++++++++++++++++++++++++++++++++++++++++++
	//  Input Tree
	//++++++++++++++++++++++++++++++++++++++++++
	TFile *f=new TFile(dirfilename.c_str());
	if(!f){
		std::cerr << "ERROR: Cant find file" << std::endl;
		return 1;
	}
	
	TTree* intree = (TTree*)f->Get("anal_tree");
	int nevent=intree->GetEntries();
	std::cerr << "NumberOfEvents : " << nevent << std::endl;
	int intree_ev_num;
	//anode
	double intree_a_hg_sum_pulse_height,intree_a_lg_sum_pulse_height;
	double intree_a_hg_sum_charge,intree_a_lg_sum_charge;
	double intree_a_hg_sum_tot;
	vector<double> intree_a_hg_pulse_height,intree_a_lg_pulse_height;
	vector<double> intree_a_hg_charge,intree_a_lg_charge;
	vector<double> intree_a_hg_mainpeak_time,intree_a_lg_mainpeak_time;
	vector<double> intree_a_hg_mainrise_time;
	vector<double> intree_a_hg_mainfall_time;
	vector<double> intree_a_hg_tot;
	//cathode
	double intree_c_hg_sum_pulse_height,intree_c_lg_sum_pulse_height;
	double intree_c_hg_sum_charge,intree_c_lg_sum_charge;
	double intree_c_hg_sum_tot;
	vector<double> intree_c_hg_pulse_height,intree_c_lg_pulse_height;
	vector<double> intree_c_hg_charge,intree_c_lg_charge;
	vector<double> intree_c_hg_mainpeak_time,intree_c_lg_mainpeak_time;
	vector<double> intree_c_hg_mainrise_time;
	vector<double> intree_c_hg_mainfall_time;
	vector<double> intree_c_hg_tot;
	//detector
	vector<double> intree_xz_x,intree_xz_z,intree_yz_y,intree_yz_z;
	double intree_ave_x,intree_ave_y,intree_ave_z;
	//minority (absolute z reconstruction)
	vector<double> intree_a_hg_minopeak_time;
	vector<double> intree_c_hg_minopeak_time;
	vector<double> intree_dt;
	vector<double> intree_abs_x,intree_abs_y,intree_abs_z;
	double intree_ave_abs_x, intree_ave_abs_y, intree_ave_abs_z;

	intree->SetBranchAddress("a_hg_sum_pulse_height",&intree_a_hg_sum_pulse_height);
	intree->SetBranchAddress("a_lg_sum_pulse_height",&intree_a_lg_sum_pulse_height);
	intree->SetBranchAddress("a_hg_sum_charge",      &intree_a_hg_sum_charge);
	intree->SetBranchAddress("a_lg_sum_charge",      &intree_a_lg_sum_charge);
	intree->SetBranchAddress("a_hg_sum_tot",         &intree_a_hg_sum_tot);
	intree->SetBranchAddress("a_hg_pulse_height",    &intree_a_hg_pulse_height);
	intree->SetBranchAddress("a_lg_pulse_height",    &intree_a_lg_pulse_height);
	intree->SetBranchAddress("a_hg_charge",          &intree_a_hg_charge);
	intree->SetBranchAddress("a_lg_charge",          &intree_a_lg_charge);
	intree->SetBranchAddress("a_hg_mainpeak_time",   &intree_a_hg_mainpeak_time);
	intree->SetBranchAddress("a_lg_mainpeak_time",   &intree_a_lg_mainpeak_time);
	intree->SetBranchAddress("a_hg_mainrise_time",   &intree_a_hg_mainrise_time);
	intree->SetBranchAddress("a_hg_mainfall_time",   &intree_a_hg_mainfall_time);
	intree->SetBranchAddress("a_hg_tot",             &intree_a_hg_tot);

	intree->SetBranchAddress("c_hg_sum_pulse_height",&intree_c_hg_sum_pulse_height);
	intree->SetBranchAddress("c_lg_sum_pulse_height",&intree_c_lg_sum_pulse_height);
	intree->SetBranchAddress("c_hg_sum_charge",      &intree_c_hg_sum_charge);
	intree->SetBranchAddress("c_lg_sum_charge",      &intree_c_lg_sum_charge);
	intree->SetBranchAddress("c_lg_sum_tot",         &intree_c_hg_sum_tot);
	intree->SetBranchAddress("c_hg_pulse_height",    &intree_c_hg_pulse_height);
	intree->SetBranchAddress("c_lg_pulse_height",    &intree_c_lg_pulse_height);
	intree->SetBranchAddress("c_hg_charge",          &intree_c_hg_charge);
	intree->SetBranchAddress("c_lg_charge",          &intree_c_lg_charge);
	intree->SetBranchAddress("c_hg_mainpeak_time",   &intree_c_hg_mainpeak_time);
	intree->SetBranchAddress("c_lg_mainpeak_time",   &intree_c_lg_mainpeak_time);
	intree->SetBranchAddress("c_hg_mainrise_time",   &intree_c_hg_mainrise_time);
	intree->SetBranchAddress("c_hg_mainfall_time",   &intree_c_hg_mainfall_time);
	intree->SetBranchAddress("c_hg_tot",             &intree_c_hg_tot);

	intree->SetBranchAddress("xz_x",                 &intree_xz_x);
	intree->SetBranchAddress("xz_z",                 &intree_xz_z);
	intree->SetBranchAddress("yz_y",                 &intree_yz_y);
	intree->SetBranchAddress("yz_z",                 &intree_yz_z);
	intree->SetBranchAddress("ave_x",                &intree_ave_x);
	intree->SetBranchAddress("ave_y",                &intree_ave_y);
	intree->SetBranchAddress("ave_z",                &intree_ave_z);

	intree->SetBranchAddress("a_hg_minopeak_time",   &intree_a_hg_minopeak_time);
	intree->SetBranchAddress("c_hg_minopeak_time",   &intree_c_hg_minopeak_time);
	intree->SetBranchAddress("dt",                   &intree_dt);
	intree->SetBranchAddress("abs_x",                &intree_abs_x);
	intree->SetBranchAddress("abs_y",                &intree_abs_y);
	intree->SetBranchAddress("abs_z",                &intree_abs_z);
	intree->SetBranchAddress("ave_abs_x",            &intree_ave_abs_x);
	intree->SetBranchAddress("ave_abs_y",            &intree_ave_abs_y);
	intree->SetBranchAddress("ave_abs_z",            &intree_ave_abs_z);

	
	//++++++++++++++++++++++++++++++++++++++++++
	//  Output
	//++++++++++++++++++++++++++++++++++++++++++
	std::string filename;
	std::string::size_type pos = dirfilename.find("/");
	while(pos != std::string::npos){
		dirfilename = dirfilename.substr(pos+1);
		pos = dirfilename.find("/");
	}
	pos = dirfilename.find(".root");
	filename = dirfilename.substr(0,pos);
	std::string outfilename = filename.substr(0, index) + "_ene.root";
	TFile *fout=new TFile(outfilename.c_str(),"recreate");
	TTree* outtree = new TTree("ene_tree");
	double cal

	//++++++++++++++++++++++++++++++++++++++++++
	//  main loop
	//++++++++++++++++++++++++++++++++++++++++++
	for(int ev=0;ev<nevent;ev++){
	//for(int ev=0;ev<1;ev++){ //to test
		// Get Event Info
		intree->GetEntry(ev);
		if(ev%1==0) std::cerr << "\rEnergy Analysis Start ... : "<< ev << "/" << nevent << std::flush;
				
	}

	//++++++++++++++++++++++++++++++++++++++++++
	//  Draw
	//++++++++++++++++++++++++++++++++++++++++++

	fout->Close();
	std::cout << "" <<std::endl;
	std::cout << "--- energy end ---" <<std::endl;

	//app.Run();

	return 1;
}





