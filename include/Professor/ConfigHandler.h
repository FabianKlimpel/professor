#ifndef __CONFIGHANDLER__H
#define __CONFIGHANDLER__H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>

using namespace std;

/**
* This class handles a config file provided by the user.
* It handles the reading of the file and stores all necessary information in the object.
*/
class ConfigHandler
{
public:
	//Constructor. If no filename is given, default parameters can be set
	ConfigHandler(string configfile);
	ConfigHandler();

	/**
	* @_summaryflag, @_outdotflag: flags for writing additional summaries
	* @_covmat: flag for writing covariance matrices
	* @_thresholdfit: threshold for the RR constrain in the fitting
	* @_thresholddata: threshold for the RR constrain in the hypercube-fitting of the data
	* @_thresholderr: threshold for the RR constrain in getter of the fitting errors
	* @_chi2mean: number of modified chi2's to store in order to state a best value regarding the modified chi2
	* @_kappa: shifting parameter if the RR constrain is applied
	* @_exponent: exponent for the distance weighting of the hypercube
	*/
	int _summaryflag = 1, _outdotflag = 0, _covmat = 1;
	double _thresholdfit = 1e-10, _thresholddata = 1e-10, _thresholderr = 1e-10, _chi2mean = 100, _kappa = 0, _exponent = 1;

private:
	//setter for the member variables
	void read_thresholdfit(string line);
	void read_thresholddata(string line);
	void read_thresholderr(string line);
	void read_chi2mean(string line);
	void read_kappa(string line);
	void read_exponent(string line);
	void read_summaryflag(string line);
	void read_outdotflag(string line);
	void read_runcombs(string line);
	void read_leave_out(string line);
	void read_rngseed(string line);
	void read_covmat(string line);
};

#endif
