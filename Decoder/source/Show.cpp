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
//USER
#include "Event.h"

#define N_CHANNEL 64

#define DEBUG 0
using namespace std;

//------------- main -------------------------------------------------
int main(int argc,char *argv[]){
  clock_t t1,t2;

  if(argc <2){
    std::cerr << "Usage:" << std::endl;
    std::cerr << "/.main [filename.root]" << std::endl;
  }

  string filename=argv[1];
  std::string::size_type index = filename.find(".root");
  if( index == std::string::npos ) { 
    std::cout << "Failure!!!" << std::endl;
    return 1;
  }

  t1=clock();
  std::cerr << "======================================" << std::endl;
  std::cerr << "Read ROOT file" << std::endl;
  std::cerr << "Visualizer for 0.1c NIuTPC" << std::endl;
  std::cerr << "Version 0.1" << std::endl;
  std::cerr << "======================================" << std::endl;
  std::cerr << "file name: " << filename << endl;

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

  TCanvas *c1=new TCanvas("c1","",800,800);
  //  TPad *pad=new TPad("pad","",0,0,1,1);
  c1->cd();
  //  pad->Draw();




  //-------- roop start ------//
  //  for(int i=0;i<nevent;i++){

    tree->GetEntry(0);

    int event_id=    event->GetEventID();
    int sampling_num=event->GetSamplingNum();
    vector< vector<double> > adc=event->GetADC();
    /*
    TH2F *h_strip=new TH2F("h_strip","h_strip",N_CHANNEL,
			   0,N_CHANNEL,sampling_num,0,sampling_num);
    h_strip->SetStats(0);
    for(int j=0;j<N_CHANNEL;j++){
	for(int k=0;k<sampling_num;k++){
	  h_strip->SetBinContent(j+1,k+1,adc[j][k]);
	}
    }


    cerr << "EventID: "<< event_id << endl;
    h_strip->SetStats(0);
    for(int j=0;j<N_CHANNEL;j++){
	for(int k=0;k<sampling_num;k++){
	  //	  cerr << adc[j][k] << " " ;
	}
    }


    //    pad->cd();
    h_strip->Draw("lego");
    */
    TH1F *h=new TH1F("h","h",sampling_num,0,sampling_num);
    for(int k=0;k<sampling_num;k++){
      h->SetBinContent(k+1,adc[32][k]);
    }

    h->Draw();

	//      }

  //----------- End of output ------------//
  std::cerr << std::endl;
  std::cerr << "======================================" << std::endl;
  t2=clock();
  std::cerr << "execution time: " << (t2-t1)/CLOCKS_PER_SEC << std::endl;
  app.Run();
  
  return 1;
}
    
