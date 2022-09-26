
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
	std::string::size_type index = dirfilename.find("_ana.root");
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
	double hg_anode_threshold    = ni_conf->hg_anode_threshold;
	double lg_anode_threshold    = ni_conf->lg_anode_threshold;
	double hg_cathode_threshold  = ni_conf->hg_cathode_threshold;
	double lg_cathode_threshold  = ni_conf->lg_cathode_threshold;
	double minority_threshold    = ni_conf->minority_threshold;
	// double minority_ROI_start    = ni_conf->minority_ROI_start;
	// double minority_ROI_end      = ni_conf->minority_ROI_end;


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
	
	TTree* intree = (TTree*)f->Get("ana_tree");
	int nevent=intree->GetEntries();
	std::cerr << "NumberOfEvents : " << nevent << std::endl;
	int intree_ev_num;
	//anode
	double intree_a_hg_sum_pulse_height,intree_a_lg_sum_pulse_height;
	double intree_a_hg_sum_charge,intree_a_lg_sum_charge;
	double intree_a_hg_sum_tot;
	vector<double> *intree_a_hg_pulse_height=0,*intree_a_lg_pulse_height=0;
	vector<double> *intree_a_hg_charge=0,*intree_a_lg_charge=0;
	vector<double> *intree_a_hg_mainpeak_time=0,*intree_a_lg_mainpeak_time=0;
	vector<double> *intree_a_hg_mainrise_time=0;
	vector<double> *intree_a_hg_mainfall_time=0;
	vector<double> *intree_a_hg_tot=0;
	//cathode
	double intree_c_hg_sum_pulse_height,intree_c_lg_sum_pulse_height;
	double intree_c_hg_sum_charge,intree_c_lg_sum_charge;
	double intree_c_hg_sum_tot;
	vector<double> *intree_c_hg_pulse_height=0,*intree_c_lg_pulse_height=0;
	vector<double> *intree_c_hg_charge=0,*intree_c_lg_charge=0;
	vector<double> *intree_c_hg_mainpeak_time=0,*intree_c_lg_mainpeak_time=0;
	vector<double> *intree_c_hg_mainrise_time=0;
	vector<double> *intree_c_hg_mainfall_time=0;
	vector<double> *intree_c_hg_tot=0;
	//detector
	vector<double> *intree_xz_x=0,*intree_xz_z=0,*intree_yz_y=0,*intree_yz_z=0;
	double intree_ave_x,intree_ave_y,intree_ave_z;
	//minority (absolute z reconstruction)
	vector<double> *intree_a_hg_minopeak_time=0;
	vector<double> *intree_c_hg_minopeak_time=0;
	vector<double> *intree_dt=0;
	vector<double> *intree_abs_x=0,*intree_abs_y=0,*intree_abs_z=0;
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
	intree->SetBranchAddress("c_hg_sum_tot",         &intree_c_hg_sum_tot);
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
	//  visualize
	//++++++++++++++++++++++++++++++++++++++++++
	//energy
	TH1D* h_sum_hg_pulse_height = new TH1D("h_sum_hg_pulse_height","h_sum_hg_pulse_height",40,0,100000);
	TH1D* h_sum_hg_charge       = new TH1D("h_sum_hg_charge"      ,"h_sum_hg_charge"      ,40,0,5000000);
	TH1D* h_sum_lg_pulse_height = new TH1D("h_sum_lg_pulse_height","h_sum_lg_pulse_height",40,0,5000);
	TH1D* h_sum_lg_charge       = new TH1D("h_sum_lg_charge"      ,"h_sum_lg_charge"      ,40,0,250000);
	//position
	TH2D* h_xy                  = new TH2D("h_xy"                 ,"h_xy"                 ,64,-1.4,1.16,64,-1.4,1.16); //bin width = 400um
	TH2D* h_xz                  = new TH2D("h_xz"                 ,"h_xz"                 ,64,-1.4,1.16,360,0,14.4); //bin width = 400um
	TH2D* h_yz                  = new TH2D("h_yz"                 ,"h_yz"                 ,360,0,14.4,64,-1.4,1.16); //bin width = 400um
	//gamma
	TH2D* h_a_sum_tot           = new TH2D("h_a_sum_tot"          ,"h_a_sum_tot"          ,40,0,100000,64,0,6000); //us
	TH2D* h_a_ave_tot           = new TH2D("h_a_ave_tot"          ,"h_a_ave_tot"          ,40,0,100000,64,0,200); //us
	//direction
	//minority
	TH1D* h_dt                  = new TH1D("h_dt"                 ,"h_dt"                 ,1000,0,400); //bin width 0.4us
	TH1D* h_abs_z               = new TH1D("h_abs_z"              ,"h_abs_z"              ,64,0,14.4); //bin width = 400um
	TH1D* h_strip_dt            = new TH1D("h_strip_dt"           ,"h_strip_dt"           ,1000,0,400); //bin width 0.4us
	TH1D* h_strip_abs_z         = new TH1D("h_strip_abs_z"        ,"h_strip_abs_z"        ,64,0,14.4); //bin width = 400um
	TH2D* h_abs_xz              = new TH2D("h_abs_xz"             ,"h_abs_xz"             ,64,-1.4,1.16,360,0,14.4); //bin width = 400um
	TH2D* h_abs_yz              = new TH2D("h_abs_yz"             ,"h_abs_yz"             ,360,0,14.4,64,-1.4,1.16); //bin width = 400um


	//++++++++++++++++++++++++++++++++++++++++++
	//  main loop
	//++++++++++++++++++++++++++++++++++++++++++
	for(int ev=0;ev<nevent;ev++){
	//for(int ev=0;ev<1;ev++){ //to test
		// Get Event Info
		intree->GetEntry(ev);
		if(ev%1==0) std::cerr << "\rAnalysis Start ... : "<< ev << "/" << nevent << std::flush;
		
		//main peak analysis
		for(int i=0;i<intree_a_hg_mainpeak_time->size();i++){
			h_xz->Fill(intree_xz_x->at(i),intree_xz_z->at(i));
			h_abs_xz->Fill(intree_xz_x->at(i),(intree_xz_z->at(i)-intree_ave_z)+intree_ave_abs_z);
		}
		for(int i=0;i<intree_c_hg_mainpeak_time->size();i++){
			h_yz->Fill(intree_yz_z->at(i),intree_yz_y->at(i));
			h_abs_yz->Fill((intree_yz_z->at(i)-intree_ave_z)+intree_ave_abs_z,intree_yz_y->at(i));
		}
		//minority peak analysis
		for(int i=0;i<intree_a_hg_minopeak_time->size();i++){
			h_strip_dt->Fill(intree_dt->at(i));
			h_strip_abs_z->Fill(intree_abs_z->at(i));
		}
		//++++++++++++++++++++++++++++++++++++++++++
		//  Fill
		//++++++++++++++++++++++++++++++++++++++++++
		//spectrum
		h_sum_hg_pulse_height->Fill(intree_a_hg_sum_pulse_height);
		h_sum_lg_pulse_height->Fill(intree_a_lg_sum_pulse_height);
		h_sum_hg_charge->Fill(intree_a_hg_sum_charge);
		h_sum_lg_charge->Fill(intree_a_lg_sum_charge);
		//tot
		h_a_sum_tot->Fill(intree_a_hg_sum_pulse_height,intree_a_hg_sum_tot);
		h_a_ave_tot->Fill(intree_a_hg_sum_pulse_height,intree_a_hg_sum_tot/intree_a_hg_tot->size());
		//position
		h_xy->Fill(intree_ave_x,intree_ave_y);
		//minority
		h_dt->Fill(std::accumulate(intree_dt->begin(),intree_dt->end(),0.0)/intree_dt->size());
		h_abs_z->Fill(intree_ave_abs_z);
	}

	//++++++++++++++++++++++++++++++++++++++++++
	//  Draw
	//++++++++++++++++++++++++++++++++++++++++++
	gStyle->SetPalette(kRainBow);
	TCanvas* c_vis = new TCanvas("c_vis","c_vis",0,0,1200,900);
	c_vis->Divide(4,3);
	//energy
	c_vis->cd(7);
	h_sum_hg_pulse_height->GetXaxis()->SetTitle("sum_pulse_height(mV)");
	h_sum_hg_pulse_height->GetYaxis()->SetTitle("counts");
	h_sum_hg_pulse_height->Draw();
	c_vis->cd(11);
	h_sum_hg_charge->GetXaxis()->SetTitle("sum_charge");
	h_sum_hg_charge->GetYaxis()->SetTitle("counts");
	h_sum_hg_charge->Draw();
	c_vis->cd(8);
	h_sum_lg_pulse_height->GetXaxis()->SetTitle("sum_pulse_height(mV)");
	h_sum_lg_pulse_height->GetYaxis()->SetTitle("counts");
	h_sum_lg_pulse_height->Draw();
	c_vis->cd(12);
	h_sum_lg_charge->GetXaxis()->SetTitle("sum_charge");
	h_sum_lg_charge->GetYaxis()->SetTitle("counts");
	h_sum_lg_charge->Draw();
	//position
	c_vis->cd(5);
	h_xz->GetXaxis()->SetTitle("x(cm)");
	h_xz->GetYaxis()->SetTitle("z(cm)");
	h_xz->Draw("colz");
	c_vis->cd(9);
	h_xy->GetXaxis()->SetTitle("x(cm)");
	h_xy->GetYaxis()->SetTitle("y(cm)");
	h_xy->Draw("colz");
	c_vis->cd(10);
	h_yz->GetXaxis()->SetTitle("z(cm)");
	h_yz->GetYaxis()->SetTitle("y(cm)");
	h_yz->Draw("colz");
	//tot
	c_vis->cd(2);
	h_a_sum_tot->GetXaxis()->SetTitle("sum_a_pulse_hight");
	h_a_sum_tot->GetYaxis()->SetTitle("sum_a_tot");
	h_a_sum_tot->Draw("colz");
	c_vis->cd(6);
	h_a_ave_tot->GetXaxis()->SetTitle("sum_a_pulse_hight");
	h_a_ave_tot->GetYaxis()->SetTitle("ave_a_tot");
	h_a_ave_tot->Draw("colz");
	

	TCanvas* c_mino = new TCanvas("c_mino","c_mino",0,0,1200,900);
	c_mino->Divide(4,3);
	c_mino->cd(9);
	h_xy->Draw("colz");
	c_mino->cd(5);
	h_abs_xz->GetXaxis()->SetTitle("x(cm)");
	h_abs_xz->GetYaxis()->SetTitle("z(cm)");
	h_abs_xz->Draw("colz");
	c_mino->cd(10);
	h_abs_yz->GetXaxis()->SetTitle("z(cm)");
	h_abs_yz->GetYaxis()->SetTitle("y(cm)");
	h_abs_yz->Draw("colz");
	c_mino->cd(7);
	h_strip_dt->GetXaxis()->SetTitle("main - minority(us)");
	h_strip_dt->GetYaxis()->SetTitle("counts");
	h_strip_dt->Draw();
	c_mino->cd(11);
	h_strip_abs_z->GetXaxis()->SetTitle("absolute z(cm)");
	h_strip_abs_z->GetYaxis()->SetTitle("counts");
	h_strip_abs_z->Draw();
	c_mino->cd(8);
	h_dt->GetXaxis()->SetTitle("main - minority(us)");
	h_dt->GetYaxis()->SetTitle("counts");
	h_dt->Draw();
	c_mino->cd(12);
	h_abs_z->GetXaxis()->SetTitle("absolute z(cm)");
	h_abs_z->GetYaxis()->SetTitle("counts");
	h_abs_z->Draw();
	


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
	std::string outfilename = filename.substr(0, index) + "_vis.root";
	TFile *fout=new TFile(outfilename.c_str(),"recreate");
	h_sum_hg_pulse_height->Write();
	h_sum_lg_pulse_height->Write();
	h_sum_hg_charge->Write();
	h_sum_lg_charge->Write();
	h_xy->Write();
	h_xz->Write();
	h_yz->Write();
	h_a_sum_tot->Write();
	h_a_ave_tot->Write();
	h_abs_xz->Write();
	h_abs_yz->Write();
	h_strip_dt->Write();
	h_strip_abs_z->Write();
	h_dt->Write();
	h_abs_z->Write();
	fout->Close();
	std::cout << "" <<std::endl;
	std::cout << "--- vis end ---" <<std::endl;

	app.Run();

	return 1;
}




