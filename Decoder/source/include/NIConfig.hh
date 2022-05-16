#ifndef NIConfig_hh
#define NIConfig_hh

#include <stdio.h>
#include <string>
#include <iostream>

class NIConfig
{
private:

public:
    NIConfig();
    ~NIConfig();

    int    offset_sampling = 500; //clock
    double cal_factor   = 0.000237; //keV/ADC
    double driftV_main  = 8.5; //cm/us
    double driftV_mino  = 8.9; //cm/us
    double tot_anode_threshold = 50; //mV
    double lg_anode_threshold = 5; //mV
    double tot_cathode_threshold = -50; //mV
    double lg_cathode_threshold = -5; //mV
    double minority_threshold = 40; //mV
    double minority_ROI_start = 50; //us
    double minority_ROI_end   = 500; //us
		
    bool ReadConfigJSON(std::string filename);
    void PrintConfigJSON();
};

#endif // NIConfig_hh
