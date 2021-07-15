//============================================================
// Decorder for NIuTPC Logger
//------------------------------------------------------------
// VERSION 
// ------------------------------------------------------------
//
// Input data format is 
//   HEADER(16byte) + ADCDATA(4000*2*64) + TIMESTAMP(4byte)
//------------------------------------------------------------
// Update : 19.Jan 2018
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
#include <algorithm>

//ROOT
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
// USER
// parameters
#define HEADER_SIZE 16 //byte
#define N_CHANNEL 64
#define N_SAMPLE 4000  // always 4000
#define SAMPLING_HELZ 2.5e6 //default
using namespace std;

int vector_finder(std::vector<int> vec, int number) {
  vector<int>::iterator itr = std::find(vec.begin(), vec.end(), number);
  size_t index = std::distance( vec.begin(), itr );
  if (index != vec.size()) {
    return 1;
  }
  else {
    return 0;
  }
}



//===========================================================================
int main(int argc,char *argv[]){
  
  if(argc != 2) {
    std::cout << "<ERROR> wrong input argument"
	      << std::endl;
    return 1;
  }

  string filename=argv[1];
  std::string::size_type index = filename.find(".dat");
  if( index == std::string::npos ) { 
    std::cout << "FILE ERROR!!!" << std::endl;
    return 1;
  }

  // ------- get file size ------//
  FILE *fp;
  long file_size;
  char *buffer;
  struct stat stbuf;
  int fd;
  fd = open(filename.c_str(), O_RDONLY);
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
  ifstream fin( filename.c_str(), ios::in | ios::binary );
  if (!fin){
    cerr << "ERROR: Cannot open file" << endl;
    return 1;
  }

  //parameters
  int event_id=0; // add event id infomation
  int current_size=0; 
  int event_num=0;
  int trigger_num=0;

  vector<int> vec_trig_0,vec_trig_1,vec_trig_2,vec_trig_3;
  int max_trig_num=0;
  //------- Start Read -------//
  while(!fin.eof()){
  // header read
    unsigned char header[HEADER_SIZE];
    fin.read( ( char * ) &header[0], HEADER_SIZE*sizeof( char ) );

    unsigned int* p = (unsigned int *)&header;
    int m_magic = p[0];

    // ID
    int module_num = ntohl(p[1]);
    module_num = module_num & 0xffffff;


    
    // DATA length
    int length = ntohl(p[2]);// unit?
    
    // trigger count
    int trigger = ntohl(p[3]);

    if(max_trig_num < trigger) max_trig_num=trigger;
    if(module_num == 0) vec_trig_0.push_back(trigger);
    if(module_num == 1) vec_trig_1.push_back(trigger);
    if(module_num == 2) vec_trig_2.push_back(trigger);
    if(module_num == 3) vec_trig_3.push_back(trigger);
    
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

    unsigned short* short_p=(unsigned short *)&data_s[0];
    
    cout << " Trigger : " << trigger << " Module     : " << module_num
	 << " Length  : " << length  << " Time Stamp : " << timestamp
	 << " 1st ADC : " << ntohs(*short_p)
	 << endl;


  }

	for(int i=0;i<vec_trig_2.size();i++){
		cout << vec_trig_2.at(i) <<endl;
	}

  int trig_num_0=vec_trig_0.size();
  int trig_num_1=vec_trig_1.size();
  int trig_num_2=vec_trig_2.size();
  int trig_num_3=vec_trig_3.size();
  int coin_num=0;
  for(int i=0;i<vec_trig_0.size();i++){
    bool coin_flag0=0;
    bool coin_flag1=0;
    bool coin_flag2=0;
		int ret_1=vector_finder(vec_trig_1,vec_trig_0.at(i));
    if( ret_1 == 1 ) coin_flag0=1; // found
    if( coin_flag0 ){
			int ret_2=vector_finder(vec_trig_2,vec_trig_0.at(i));
    	if( ret_2 ==1 ) coin_flag1=1;
		}
    if( coin_flag1 ){
			int ret_3=vector_finder(vec_trig_3,vec_trig_0.at(i));
  		if( ret_3 == 1 ) coin_flag2 = 1;
		}
		if(coin_flag2)coin_num++;
	}

  cout << "Total trigger number : " << max_trig_num << endl;
  cout << "Module 0 trigger number  : " << trig_num_0 << endl;
  cout << "Module 1 trigger number  : " << trig_num_1 << endl;    
  cout << "Module 2 trigger number  : " << trig_num_2 << endl;
  cout << "Module 3 trigger number  : " << trig_num_3 << endl;    
  cout << "Coincidence number  : " << coin_num << endl;      
  cout << "Efficiency  : " << (double)coin_num/((double)max_trig_num+1) << endl;      
  cerr << endl;


  fin.close();
  
  return 0;


}
