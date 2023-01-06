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
	printf("ANODE HG Threshold               : %.1lf mV\n",hg_anode_threshold);
	printf("ANODE LG Threshold               : %.1lf mV\n",lg_anode_threshold);
	printf("CATHODE HG Threshold             : %.1lf mV\n",hg_cathode_threshold);
	printf("CATHODE LG Threshold             : %.1lf mV\n",lg_cathode_threshold);
	printf("Minority Threshold               : %.1lf mV\n",minority_threshold);
	printf("Minority ROI range               : %.1lf us\n",minority_ROI_range);
	printf("Minority ROI offset              : %.1lf us\n",minority_ROI_offset);
	printf("Is Alpha Calibration?            : %d (1: yes, 0: no)\n",is_alpha_calib);
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

	if(boost::optional<double> buf = pt.get_optional<double>("config.hg_anode_threshold")){
		this->hg_anode_threshold = pt.get<double>("config.hg_anode_threshold");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.lg_anode_threshold")){
		this->lg_anode_threshold = pt.get<double>("config.lg_anode_threshold");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.hg_cathode_threshold")){
		this->hg_cathode_threshold = pt.get<double>("config.hg_cathode_threshold");
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

	if(boost::optional<double> buf = pt.get_optional<double>("config.minority_ROI_range")){
		this->minority_ROI_range = pt.get<double>("config.minority_ROI_range");
	}else{
	}

	if(boost::optional<double> buf = pt.get_optional<double>("config.minority_ROI_offset")){
		this->minority_ROI_offset = pt.get<double>("config.minority_ROI_offset");
	}else{
	}

	if(boost::optional<int> buf = pt.get_optional<int>("config.is_alpha_calib")){
		this->is_alpha_calib = pt.get<int>("config.is_alpha_calib");
	}else{
	}

    return true;
}
