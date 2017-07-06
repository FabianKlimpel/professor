#ifndef __GRADHANDLER__H
#define __GRADHANDLER__H

#include <iostream>
#include <stdlib.h>
#include <vector>
#include "Professor/LinAlg.h"
#include "Professor/ParamPoints.h"

//forward declaration
class ParamPoints;

using namespace std;
/**
* This class handles the gradient calculation of the anchor points
*/
class GradHandler
{
//allowing FitHandler to access @_gradvecs
friend class FitHandler;

public:

	//Constructor
	GradHandler();
	
	//calculates the gradient vectors for the smoothness convergence check
	void calculategradvectors(Professor::ParamPoints& pts, const vector<double>& ptvals, double& threshold, double& power, double& kappa);

	//calculates the dot product between every pair of gradient vectors
	vector<double> getallgraddotproducts();

private:
	
	/**
	* @_gradvecs: storage of the gradient vectors
	*/
	vector<vector<double>> _gradvecs;

	//calculates the normalvector for an anchor point
	vector<double> getgrad(Professor::ParamPoints& pts, const vector<double>& ptvals, size_t i, double threshold, int power, double kappa);

	//calculates the fit parameters
	vector<double> getfitparams(vector<vector<double> >& points, vector<double>& bpoints, double threshold, int power, double kappa);

	//creates a QR-decomposition
	void expandqrlocal(vector<vector<double> >& mlocal, vector<vector<double> >& qlocal, vector<vector<double> >& rlocal, size_t i);
};

#endif

