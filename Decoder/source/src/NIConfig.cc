#include "NIConfig.hh"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>


NIConfig::NIConfig()
{
}

NIConfig::~NIConfig()
{
}


void NIConfig::PrintConfigJSON()
{
	
	printf("---------- configuration parameters ----------\n");
	printf("Waveform Offset Calc Sampling    : %d clock\n",offset_sampling);
	printf("Calibration Factor               : %lf keV/ADC\n",cal_factor);
	printf("Drift Velocity (Main Charge)     : %.2lf cm/us\n",driftV_main);
	printf("Drift Velocity (Minority Charge) : %.2lf cm/us\n",driftV_mino);
	printf("ANODE TOT(HG) Threshold          : %.1lf mV\n",tot_anode_threshold);
	printf("ANODE LG Threshold               : %.1lf mV\n",lg_anode_threshold);
	printf("CATHODE TOT(HG) Threshold        : %.1lf mV\n",tot_cathode_threshold);
	printf("CATHODE LG Threshold             : %.1lf mV\n",lg_cathode_threshold);
	printf("Minority Threshold               : %.1lf mV\n",minority_threshold);
	printf("Minority Saerch ROI              : %.1lf ~ %.1lf us\n",minority_ROI_start,minority_ROI_end);
	printf("----------------------------------------------\n");

}

bool NIConfig::ReadConfigJSON(std::string conffilename)
{
	std::string::size_type index_conf = conffilename.find(".json");
	if( index_conf == std::string::npos ) { 
		std::cout << "Failure!!!" << std::endl;
        return false;
	}

	boost::property_tree::ptree pt;
	read_json(conffilename,pt);

	if(boost::optional<int> buf = pt.get_optional<int>("config.offset_sampling")){
		this->offset_sampling = pt.get<int>("config.offset_sampling");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.cal_factor")){
		this->cal_factor = pt.get<double>("config.cal_factor");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.driftV_main")){
		this->driftV_main = pt.get<double>("config.driftV_main");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.driftV_mino")){
		this->driftV_mino = pt.get<double>("config.driftV_mino");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.tot_anode_threshold")){
		this->tot_anode_threshold = pt.get<double>("config.tot_anode_threshold");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.lg_anode_threshold")){
		this->lg_anode_threshold = pt.get<double>("config.lg_anode_threshold");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.tot_cathode_threshold")){
		this->tot_cathode_threshold = pt.get<double>("config.tot_cathode_threshold");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.lg_cathode_threshold")){
		this->lg_cathode_threshold = pt.get<double>("config.lg_cathode_threshold");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.minority_threshold")){
		this->minority_threshold = pt.get<double>("config.minority_threshold");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.minority_ROI_start")){
		this->minority_ROI_start = pt.get<double>("config.minority_ROI_start");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.minority_ROI_end")){
		this->minority_ROI_end = pt.get<double>("config.minority_ROI_end");
	}else{
	}

    return true;
}
