#ifndef __FITHANDLER__H
#define __FITHANDLER__H

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <math.h>
#include "Professor/LinAlg.h"
#include <omp.h>
#include "Professor/OutputHandler.h"
#include <limits>
#include "Professor/ParamPoints.h"
#include "Professor/GradHandler.h"
#include "Professor/Ipol.h"
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//forward declaration
class OutputHandler;

/**
* This class handles the fitting procedure.
*/
class FitHandler
{

public:
	//Default Constructor for fitting
	FitHandler();

	//Constructor used for error calculation
	FitHandler(vector<double>& fitparams, const vector<double>& pterrs, Professor::ParamPoints& pts, int& order, int& num_ipol);

	//calculates the next iterationstep
	void nextstep(Professor::ParamPoints& pts, double& threshold, double& kappa);

	//calculator of fit errors
	void setfiterrors(size_t num_ipol);

	//getter of the smoothness
	const double getDsmooth(Professor::ParamPoints& pts, GradHandler& gh);

	//getter of Chi2
	const double getchi2();

	//getter of @_iterationcounter
	size_t getiterationcounter() {return _iterationcounter;}

	//setup of a new bin
	void startbin(Professor::ParamPoints& pts, const vector<double>& ptvals, const vector<double>& pterrs, int& num_ipol, double& threshold, double& kappa);

	//getter of the number of fitparameters
	size_t getnumfitparams();

	//getter of the fitparameters
	vector<double> getfitparams() {return _bfp;}

	//getter of @_a
	const vector<double> geta() {return _a;}

	//setter of a specific bin in a specific iterationstep
	void setiteration(Professor::ParamPoints& pts, const vector<double>& ptvals, const vector<double>& pterrs, int num_ipol, double threshold, double kappa, int bestiteration);

	//getter of the mean relative error; not in use right now
	//const double getmeanerror();

	//calculates the dot product between every gradient vector
	vector<double> getallgraddotproducts(Professor::ParamPoints& pts, GradHandler& gh);

	//calculates the error at a given point; not in use right now
	//const double geterrorat(int dim);

	//calculates the reduced Chi^2
	const double getchi2red();

	//adds 0's as fit parameters and changes the order of the fit parameters in order to become usable in Professor 2.2.1
	void sortfitparams(vector<int>& match_index);

	//getter of @_max
	int getmax() {return _max;}

	//getter of the fit errors
	vector<double> getfiterrors() {return _bfperr;}

	//getter of @_power
	vector<vector<int>> getpowers() {return _power;}

private:
	/**
	* @_m: storage of the matrix M
	* @_q, @r: storage of the matrices Q & R of the QR-decomposition
	* @_d: identity matrix
	* @_mprime: coloumnwise normalized M
	* @_b: reference values
	* @_a: vector for the fitparameters; is used as indicator of RR constrained parameters
	* @_bprime: normalized @_b
	* @_bfp, @_bfperr: best fit parameters and corresponding errors
	* @_sigma: error of the reference values
	* @_max: maximum of powers
	* @_iterationcounter: number of iterations in the fitting process
	* @_power: list that contains the powers of the variables at the terms of the fitting function
	* @_dsmooth_current: Storage of the dsmooth value in a certain step in order to avoid recalculation. The default value indicates an unset value.
	*/
	vector<vector<double>> _m, _r, _q, _d, _mprime;
	vector<double> _b, _a, _bprime, _bfp, _bfperr, _sigma;
	int _max; 
	size_t _iterationcounter;
	vector<vector<int>> _power;
	double _dsmooth_current = 2;
	

	//adding components to @_q & @_r
	void expandqr();

	//check for RR constraints
	void collinearity(double threshold, size_t i, double kappa);

	//calculate @_mprime from @_m
	void makemprime();

	//adding elements to @_m
	void increasem(Professor::ParamPoints& pts);

	//rescaling @_bfp because of the normalization
	void rescalebestfitparameters();
	
	//calculates the gradients of the function
	vector<double> gradvecfunction(Professor::ParamPoints& pts, size_t i);
	
	//this function walks through an iteration without any checks
	void nextstep_walkthrough(Professor::ParamPoints& pts);	

	//getter of a single element of the precision matrix
	double getinvcovmatelement(size_t i, size_t j);
	
};

#endif

