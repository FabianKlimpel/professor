#include "Professor/ConfigHandler.h"

/**
 * Constructor that reads a config file and sets all member variables according to the setted ones
 * @configfile: name of the config file
 * @ifile: input stream for the config file
 * @line: string that gets the content of a line in the configfile
 */
ConfigHandler::ConfigHandler(string configfile){
	
	cout << "reading config file ...";
	ifstream ifile;
	ifile.open(configfile);
	string line;
	
	if(ifile.is_open())
	{
		//linewise file reading
		while(getline(ifile, line))
		{
			//checking for signal words and calling the respective function
			if(line.substr(0, 12) == "thresholdfit")
				read_thresholdfit(line);
				
			if(line.substr(0, 13) == "thresholddata")
				read_thresholddata(line);
				
			if(line.substr(0, 12) == "thresholderr")
				read_thresholderr(line);
				
			if(line.substr(0, 8) == "chi2mean")
				read_chi2mean(line);
				
			if(line.substr(0, 5) == "kappa")
				read_kappa(line);
				
			if(line.substr(0, 8) == "exponent")
				read_exponent(line);
				
			if(line.substr(0, 7) == "summary")
				read_summaryflag(line);
				
			if(line.substr(0, 6) == "outdot")
				read_outdotflag(line);
				
			if(line.substr(0, 6) == "covmat")
				read_covmat(line);
		}		
	}

	cout << "complete" << endl;
	
}

/**
 * Default constructor. If no config file is delivered, a default configuration will be set.
 */
ConfigHandler::ConfigHandler(){
	
	_summaryflag = 1;
	_outdotflag = 0;
	_covmat = 1;
	_thresholdfit = 1e-10;
	_thresholddata = 1e-10;
	_thresholderr = 1e-10;
	_kappa = 1e-10;
	_chi2mean = 100;
	_exponent = 1;
	
}

/**
 * This function reads the threshold for the RR constrain of the fit
 * @line: contains the threshold value
 */
void ConfigHandler::read_thresholdfit(string line){
		_thresholdfit = atof(line.substr(13).c_str());
}

/**
 * This function reads the threshold for the RR constrain of the hypercube fitting
 * @line: contains the threshold value
 */
void ConfigHandler::read_thresholddata(string line){
		_thresholddata = atof(line.substr(14).c_str());	
}

/**
 * This function reads the threshold for the RR constrain of the error of the fit
 * @line: contains the threshold value
 */
void ConfigHandler::read_thresholderr(string line){
		_thresholderr = atof(line.substr(13).c_str());	
}

/**
 * This function reads the number of Chi2 values for the mean
 * @line: contains the number of points
 */
void ConfigHandler::read_chi2mean(string line){
		_chi2mean = atof(line.substr(9).c_str());	
}

/**
 * This function reads the kappa for the RR constrain
 * @line: contains the value
 */
void ConfigHandler::read_kappa(string line){
		_kappa = atof(line.substr(6).c_str());
}

/**
 * This function reads the exponent for the distance weighting of the hypercube fitting
 * @line: contains the exponent
 */
void ConfigHandler::read_exponent(string line){
		_exponent= atof(line.substr(9).c_str());	
}

/**
 * This function reads the flag for writing the summary output
 * @line: contains the value
 */
void ConfigHandler::read_summaryflag(string line){
	_summaryflag = atoi(line.substr(8).c_str());	
}

/**
 * This function reads the flag for writing the dot product summary
 * @line: contains the value
 */
void ConfigHandler::read_outdotflag(string line){
	_outdotflag = atoi(line.substr(7).c_str());
}

/**
 * This function reads the flag for writing the covariance matrices
 * @line: string that contains the flag
 */
void ConfigHandler::read_covmat(string line){
	_covmat = atoi(line.substr(7).c_str());
}
