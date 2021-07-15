//==================================================================
// Event Class for NIuTPC Data
//
// DAQ version is mujirushi(NIuTPC)
//------------------------------------------------------------------
// Author : T.Ikeda
// Update : 10. July 2016
//==================================================================

// ROOT 
#include <TObject.h>
//#include <RootInt.h>

#include <vector>
using namespace std;

class Event: public TObject
{
 private:

 //-------- HEADERS ---------//
  int module_num;
  int trigger;
  //
  int sampling_num;
  double sampling_helz;

 //-------  ADC DATA ---------//
  vector< vector<double> > a_hg_adc;
  vector< vector<double> > a_lg_adc;
  vector< vector<double> > c_hg_adc;
  vector< vector<double> > c_lg_adc;

 //------ TIME STAMP--------//
  int timestamp;

  //Event ID
 public:
  Event();
  ~Event();
 //-------- set fuctions --------//
  //void SetHeaders(int trigger,int timestamp,int sampling_num,double sampling_helz);
  void SetHeaders(int trigger,int timestamp);
  void SetAnodeHGADC(vector< vector<double> > a_hg_adc);
  void SetAnodeLGADC(vector< vector<double> > a_lg_adc);
  void SetCathodeHGADC(vector< vector<double> > c_hg_adc);
  void SetCathodeLGADC(vector< vector<double> > c_lg_adc);

 //-------- get fuctions --------//
  int GetTrigger();
  int GetTimeStamp();
  int GetSamplingNum();
  double GetSamplingHelz();

	vector<vector<double>> GetAnodeHGADC();	
	vector<vector<double>> GetAnodeLGADC();	
	vector<vector<double>> GetCathodeHGADC();	
	vector<vector<double>> GetCathodeLGADC();	

  ClassDef(Event,1);  
};
