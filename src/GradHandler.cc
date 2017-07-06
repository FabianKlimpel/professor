#include "Professor/GradHandler.h"

/**
 * Constructor
 */
GradHandler::GradHandler(){}

/**
 * This function calculates the gradient vectors of every reference data point
 * @pts: container of the anchor points
 * @ptvals: bin value of every anchor point
 * @threshold: RR constraint parameter for the fitting in the hypercubes
 * @power: specifies the power of the distance weighting
 * @kappa: shift in case of applied RR constraint
 */
void GradHandler::calculategradvectors(Professor::ParamPoints& pts, const vector<double>& ptvals, double& threshold, double& power, double& kappa){

	_gradvecs.resize(pts.numPoints());
		
	//walk over every point
	for(size_t i = 0; i < pts.numPoints(); i++)
		//calculate its gradient vector
		_gradvecs[i] = getgrad(pts, ptvals, i, threshold, power, kappa);	
}

/**
 * This function calculates the gradient vector of one point
 * @pts: container of the anchor points
 * @ptvals: bin value of every anchor point
 * @i: index of the point of interest
 * @threshold: RR constraint threshold for the fitting
 * @power: specifies the power of the distance weighting
 * @kappa: shift in case of applied RR constraint
 * @fitparams: vector that will store the fitting result and after that the gradient
 * @bpoints: list of simulated reference data that belong to the hypercube
 * @points: list of parameter values that belong to the hypercube
 * @length: storage for the length of a vector
 */
vector<double> GradHandler::getgrad(Professor::ParamPoints& pts, const vector<double>& ptvals, size_t i, double threshold, int power, double kappa){

	vector<double> fitparams, bpoints;
	fitparams.resize(pts.dim() + 1);
	vector<vector<double>> points;
	
	//adding the point of interest to the lists
	points.push_back(pts.point_scaled(i));
	bpoints.push_back(ptvals[i]);
		
	//adding the points of the hypercube to the lists
	for(size_t j = 0; j < pts.gethypercube(i).size(); j++)
		if(pts.gethypercube(i)[j] < pts.numPoints())
		{
			points.push_back(pts.point_scaled(pts.gethypercube(i)[j]));
			bpoints.push_back(ptvals[pts.gethypercube(i)[j]]);
		}

	//calculate the fit parameters
	fitparams = getfitparams(points, bpoints, threshold, power, kappa);

	//calculating the gradient = the fit parameters of the 1. order only
	fitparams.erase(fitparams.begin());
	double length = LinAlg::getabs(fitparams);
	
	//if the length of the gradient vector is != 0, it will be normalized
	//if it is = 0, then the normalized gradient vector is the vector itself
	if(length != 0)
		for(size_t j = 0; j < fitparams.size(); j++)
			fitparams[j] /= length;
		
	return fitparams;
}

/**
 * This function calculates the best fit parameters for the hypercubes
 * @points: matrix that contains all points
 * @bpoints: vector with the functionvalues
 * @threshold: threshold for artificial fitparameter setting
 * @power: list containing the powers for every value at every monomial
 * @kappa: parameter shift in the case of RR constraint
 * @result: container for the best fit parameters
 * @mlocal, @qlocal, @rlocal: local equivalent to @m, @q, @r respectively
 * @distances: vector for weighting the hypercubepoints with its distance to the center
 * @sum: storage of the sum of the distances of the hypercube points to the center point
 * @weights: list of weights of the anchor points that are later used to calculate the fit function
 * @riihat, @bihat: regularization as in @collinearity()
 * @dlocal: identity matrix
 */
vector<double> GradHandler::getfitparams(vector<vector<double>>& points, vector<double>& bpoints, double threshold, int power, double kappa){

	//setting up the parametervector, setting nans for finding paramaters to fit
	vector<double> result;
	result.resize(points[0].size() + 1);
	//setting the parameters to nan = fit parameter is not calculated yet.
	result.assign(result.size(), nan("1"));
	vector<vector<double>> mlocal, qlocal, rlocal;
	mlocal.resize(points.size());	
	
	//calculating the distances of the points of the hypercube to the center
	vector<double> distances;
	distances.resize(points.size());
	
	//the zeroth component is set to 1 and therefore it won't be weighted otherwise
	distances[0] = 1.;
	double sum = 0;
	//calculate every distance and increase its "importance decrease by distance" by the power of @power
	for(size_t i = 1; i < distances.size(); i++)
	{
		distances[i] = pow(LinAlg::getdistanceofvectors(points[0], points[i]), power);
		sum += distances[i];
	}
		
	//calculate the weights
	vector<double> weights;
	
	//setting the point itself to 1, every other point is weighted based on its distance to the central point in an [0,1]-interval. That way, the central point will get the highest weight.
	weights.push_back(1);
	for(size_t i = 1; i < distances.size(); i++)
		weights.push_back(1 - distances[i] / sum);
	
	//prepare the matrices and vector to calculate a fit function
	for(size_t i = 0; i < bpoints.size(); i++)
		bpoints[i] *= weights[i];
	for(size_t i = 0; i < points.size(); i++)
		mlocal[i].push_back(weights[i]);		
	expandqrlocal(mlocal, qlocal, rlocal, 0);
	for(size_t i = 0; i < points[0].size(); i++)
	{
		for(size_t j = 0; j < points.size(); j++)
			mlocal[j].push_back(points[j][i] * weights[j]);	
		expandqrlocal(mlocal, qlocal, rlocal, i);
	}
	
	//regularize the QR-decomposition as in @collinearity()
	double riihat, bihat;
	vector<vector<double> > dlocal;
	dlocal.resize(rlocal.size());
	for(size_t i = 0; i < dlocal.size(); i++)
		dlocal[i].assign(rlocal.size(), 0.);
	for(size_t i = 0; i < dlocal.size(); i++)
		dlocal[i][i] = 1.;
			
	for(size_t i = 0; i < rlocal.size(); i++)
		if(rlocal[i][i] < threshold)
		{
			riihat = sqrt(rlocal[i][i] * rlocal[i][i] + kappa * 1.);
			bihat = (rlocal[i][i] / riihat) * bpoints[i];
			result[i] = bihat / riihat;
		}
	
	//return the fit
	return LinAlg::getbestfitparameters(rlocal, result, LinAlg::multmatvec(LinAlg::transpose(qlocal), bpoints));
}

/**
 * This function does the same as @FitHandler::expandqr() but takes the matrices as arguments and is therefore more functionally
 * @mlocal, @qlocal, @rlocal: analog to the FitHandler member variables @_m, @_q, @_r respectively
 * @i: counter of the respective iteration
 */
void GradHandler::expandqrlocal(vector<vector<double> >& mlocal, vector<vector<double> >& qlocal, vector<vector<double> >& rlocal, size_t i){
	if(qlocal.empty() || rlocal.empty())
	{
	
		qlocal.resize(mlocal.size());
		rlocal.resize(mlocal[0].size());
			
		for(size_t k = 0; k < qlocal.size(); k++)
		{
			qlocal[k].resize(mlocal[0].size());
		}
		
		for(size_t k = 0; k < rlocal.size(); k++)
		{
			rlocal[k].resize(mlocal[0].size());
		}
		
		for(size_t k = 0; k < qlocal.size(); k++)
			qlocal[k][0] = LinAlg::getcol(mlocal,0)[k] / LinAlg::getabs(LinAlg::getcol(mlocal,0));

		rlocal[0][0] = LinAlg::getabs(LinAlg::getcol(mlocal,0));
		return;
	}


	qlocal.resize(mlocal.size());
	rlocal.resize(mlocal[0].size());


	for(size_t k = 0; k < qlocal.size(); k++)
	{
	qlocal[k].resize(mlocal[0].size());
	}
	
	for(size_t k = 0; k < rlocal.size(); k++)
	{	
	rlocal[k].resize(mlocal[0].size());
	}
	

	double tmp;
	for(size_t l = 0; l <= i; l++)
	{
		tmp = 0;
		for(size_t k = 0 ; k < qlocal.size(); k++)
		{
			tmp += qlocal[k][l] * mlocal[k][i + 1];

		}
		rlocal[l][i + 1] = tmp;
	}


	vector<double> sum;
	sum.assign(qlocal.size(), 0);
	for(size_t l = 0; l < sum.size(); l++)
		for(size_t k = 0; k <= i; k++)
			sum[l] += rlocal[k][i + 1] * qlocal[l][k];


	for(size_t k = 0; k < qlocal.size(); k++)
		qlocal[k][i + 1] = mlocal[k][i + 1] - sum[k];
		
	

	rlocal[i + 1][i + 1] = LinAlg::getabs(LinAlg::getcol(qlocal,i + 1));
	
	tmp = LinAlg::getabs(LinAlg::getcol(qlocal, i + 1));
	if(tmp != 0)
		for(size_t k = 0; k < qlocal.size(); k++)
			qlocal[k][i + 1] /= tmp;
			
}

/**
 * This function calculates every dot product of every anchor point with every anchor point
 * @result: summary of all dot products
 */
vector<double> GradHandler::getallgraddotproducts(){
	
	vector<double> result;
	
	//calculate the dot product of every vector with every other
	for(size_t i = 0; i < _gradvecs.size(); i++)
		for(size_t j = 0; j < _gradvecs.size(); j++)
			//if they are the same, the calculation will be skipped
			if(i != j)
				result.push_back(LinAlg::dotproduct(_gradvecs[i], _gradvecs[j]));
	
	return result;
}
