//==================================================================
// Event Class for NIuTPC Data
//
// DAQ version is mujirushi(NIuTPC)
//------------------------------------------------------------------
// Author : T.Ikeda
// Update : 10. July 2016
//==================================================================

#include "Event.h"
#include <numeric>

ClassImp(Event);  

Event::Event()
{
  trigger = 0;
  timestamp = 0;
	sampling_num = 0;
	sampling_helz = 0;

}


Event::~Event()
{
}

//------------------------------------------------------------------
// set functions
//------------------------------------------------------------------
void Event::SetHeaders(int trigger,int timestamp)
		       //int sampling_num,double sampling_helz)
{
   this->trigger = trigger;
   this->timestamp = timestamp;
   //this->sampling_num = sampling_num;
   //this->sampling_helz = sampling_helz;
}
/*
void Event::SetADC(vector< vector<double> > adc){
  this->adc = adc;
}
*/
/*
void Event::SetEventID(int event_id){
  this->event_id=event_id;
}
*/

void Event::SetAnodeHGADC(vector<vector<double>> a_hg_adc){
	this->a_hg_adc = a_hg_adc;
}
void Event::SetAnodeLGADC(vector<vector<double>> a_lg_adc){
	this->a_lg_adc = a_lg_adc;
}
void Event::SetCathodeHGADC(vector<vector<double>> c_hg_adc){
	this->c_hg_adc = c_hg_adc;
}
void Event::SetCathodeLGADC(vector<vector<double>> c_lg_adc){
	this->c_lg_adc = c_lg_adc;
}


//------------------------------------------------------------------
// get functions
//------------------------------------------------------------------ 
int Event::GetTrigger(){
  return trigger;
}
int Event::GetTimeStamp(){
  return timestamp;
}
int Event::GetSamplingNum(){
  return sampling_num;
}
double Event::GetSamplingHelz(){
  return sampling_helz;
}
vector<vector<double>> Event::GetAnodeHGADC(){
	return a_hg_adc;
}
vector<vector<double>> Event::GetAnodeLGADC(){
	return a_lg_adc;
}
vector<vector<double>> Event::GetCathodeHGADC(){
	return c_hg_adc;
}
vector<vector<double>> Event::GetCathodeLGADC(){
	return c_lg_adc;
}



