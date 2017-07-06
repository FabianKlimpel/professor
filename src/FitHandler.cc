#include "Professor/FitHandler.h"

/**
 * Constructor
 */
FitHandler::FitHandler(){}

/**
 * Main constructor.
 * @fitparams: vector containing the fitparameters
 * @pterrs: vector containing the uncertainties of the data
 * @pts: Object containing everything about the anchor points
 * @order: order of the polynomial function
 * @num_ipol: bin number
 */
FitHandler::FitHandler(vector<double>& fitparams, const vector<double>& pterrs, Professor::ParamPoints& pts, int& order, int& num_ipol){

	//local storage of parameters
	_bfp = fitparams;
	_sigma = pterrs;
	_power = pts.getpower(order);
	_max = order;
	
	//set the size of @_m, if not done yet
	if(_m.empty())
		_m.resize(pts.numPoints());
		
	//set up @_m
	while(_m[0].size() < getnumfitparams())
		increasem(pts);
		
	//calculate the fit errors
	setfiterrors(num_ipol);
}

/**
 * This function builds the QR-decomposition or increase the matrices, if they already exist
 * @i: represents the iterationstep
 * @_q, @_r: the respective matrices of the QR-decomposition
 * @_mprime: coloumnwise rescaled @_m, c.f. @FitHandler::makemprime() 
 * @tmp: temporary storage for the increase of @_r
 * @sum: temporary storage for the increase of @_q
 */
void FitHandler::expandqr(){
	//in the first iteration, @q & @r need to be resized
	if(_q.empty() || _r.empty())
	{
		_q.resize(_mprime.size());
		_r.resize(_mprime[0].size());
		
		for(size_t k = 0; k < _q.size(); k++)
			_q[k].resize(_mprime[0].size());
		
		for(size_t k = 0; k < _r.size(); k++)
			_r[k].resize(_mprime[0].size());
		
		//calculating the first elements by using http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5342430&tag=1 on top of eq. 5
		for(size_t k = 0; k < _q.size(); k++)
			_q[k][0] = LinAlg::getcol(_mprime,0)[k] / LinAlg::getabs(LinAlg::getcol(_mprime,0));

		_r[0][0] = LinAlg::getabs(LinAlg::getcol(_mprime,0));
		return;
	}

	//in a later iteration, the size should be adapted according to the current size of @_mprime 
	_q.resize(_mprime.size());
	_r.resize(_mprime[0].size());


	for(size_t k = 0; k < _q.size(); k++)
		_q[k].resize(_mprime[0].size());
	
	
	for(size_t k = 0; k < _r.size(); k++)
		_r[k].resize(_mprime[0].size());
	
	
	//the following calculation add the new elements to @_q & @_r
	double tmp;
	for(size_t l = 0; l <= _iterationcounter - 1; l++)
	{
		tmp = 0;
		for(size_t k = 0 ; k < _q.size(); k++)
			tmp += _q[k][l] * _mprime[k][_iterationcounter];
			
		_r[l][_iterationcounter] = tmp;
	}


	vector<double> sum;
	sum.assign(_q.size(), 0);
	for(size_t l = 0; l < sum.size(); l++)
		for(size_t k = 0; k <= _iterationcounter - 1; k++)
			sum[l] += _r[k][_iterationcounter] * _q[l][k];
			

	for(size_t k = 0; k < _q.size(); k++)
		_q[k][_iterationcounter] = _mprime[k][_iterationcounter] - sum[k];


	_r[_iterationcounter][_iterationcounter] = LinAlg::getabs(LinAlg::getcol(_q, _iterationcounter));

	tmp = LinAlg::getabs(LinAlg::getcol(_q, _iterationcounter));
	for(size_t k = 0; k < _q.size(); k++)
	{
		_q[k][_iterationcounter] /= tmp;
		
		//If the absolut value of @q[k][i + 1] is 0, the result would become nan. Setting it instead to 0 "stabilizes" the calculations.
		if(std::isnan(_q[k][_iterationcounter])) 
			_q[k][_iterationcounter] = 0;
	}
}

/**
 * This function checks for the collinearity of a coloumn of the matrix @_m in comparison with the previously coloumns.
 * This leads to a small values at the corresponding diagonalterm in @_r. The "smallness" is checked and, if it's true, the according fitparameter will be artificial set in order to
 * keep the solution from diverging.
 * @threshold: this value represents the threshold, below which a diagonalterm of @r is considered as small
 * @i: this values represents the bin that is concerned
 * @riihat, @bihat: helpervariables for setting of the fitparameter in @_a
 */
void FitHandler::collinearity(double threshold, size_t i, double kappa){

	//check if the value in @r is small
	if(fabs(_r[i][i]) < threshold)
	{
		//setting the helpervariables for the creation of a smaller @a than one would get if the fit would run normally.
		//the shift has it's source in @kappa and the according term of @_d.
		double riihat = sqrt(_r[i][i] * _r[i][i] + kappa * _d[i][i]);
		double bihat = (_r[i][i] / riihat) * _b[i];
		_a[i] = bihat / riihat;
		
	}
}

/**
 * This function normalizes every coloumn of @_m by using the absolut value of the respective coloumn.
 * @coloumn: represents the coloumn of interest
 * @abs: absolut value of the coloumn
 */
void FitHandler::makemprime(){

	//setting @_mprime to the same size as @_m
	_mprime.resize(_m.size());
	for(size_t i = 0; i < _mprime.size(); i++)
		_mprime[i].resize(_m[i].size());

	//setting the absolut value of the coloumn, because while assigning it row wise, it changes after every assigned value
	double abs = LinAlg::getabs(LinAlg::getcol(_m, _iterationcounter));
	
	//setting the normalization
	for(size_t i = 0; i < _m.size(); i++)
		_mprime[i][_iterationcounter] = _m[i][_iterationcounter] / abs;

}

/**
 * This function adds a new coloumn to @_m
 * @pts: container for the anchor points
 * @size: old number of coloumns in @_m
 * @tmp: temporary storage of the product of the different powers of the values
 */
void FitHandler::increasem(Professor::ParamPoints& pts){

	//set the size of @_m, if not done yet
	if(_m.empty())
		_m.resize(pts.numPoints());
	
	//adding a new coloumn to @_m
	size_t size = _m[0].size();
	for(size_t i = 0; i < _m.size(); i++)
		_m[i].resize(size + 1);
	
	double tmp = 1;
	//walking over every row of @_m
	for(size_t j = 0; j < _m.size(); j++)
	{	
		//in every row the product of the values and the respective powers is calculated
		for(size_t i = 0; i < pts.dim(); i++)
			tmp *= pow(pts.point_scaled(j)[i], _power[_m[0].size() - 1][i]);
		
		//the new terms will be added to the last coloumn in every row and @tmp is resetted
		_m[j][size] = tmp;
		tmp = 1;
	}
}

/**
 * This function rescales the best fit parameters due to normalization of @_m and @_b
 * @vec: vector, containing the best fit parameters
 */
void FitHandler::rescalebestfitparameters(){

	//backwards calculating of the norms that influenzed @_m and @_b in order to get the regular best fit parameters
	for(size_t i = 0; i < _bfp.size(); i++)
		_bfp[i] *= LinAlg::getabs(_b) / LinAlg::getabs(LinAlg::getcol(_m, i));
}

/**
 * This function calculates the Chi^2 value of the fit.
 * @result: resulting Chi^2
 * @functionvalue: functionvalue of the fit at a certain point
 */
const double FitHandler::getchi2(){

	double result = 0, functionvalue = 0;

	//walk over every polynomial
	for(size_t i = 0; i < _m.size(); i++)
	{
		//if the bin entry is 0, it can be a simple misscalculation in Rivet etc.
		if(_b[i] == 0)
			continue;
		
		for(size_t j = 0; j < _m[i].size(); j++)
			functionvalue += _m[i][j] * _bfp[j];
			
		//get the difference between the functionvalue and the datapoint		
		functionvalue -= _b[i];
		
		//if the uncertainty is 0, the calculation would rise a nan
		//this is prevented by setting it to the arbitrary value of 1e-10
		if(_sigma[i] == 0)
			_sigma[i] = 1e-10;
		
		//calculating the squares of @functionvalue and @sigma, divide them and add this new summand to new overall Chi^2 value @result
		result += (functionvalue * functionvalue) / (_sigma[i] * _sigma[i]);

		functionvalue = 0;
	}
	return result;
}

/**
 * This function calculates the Chi^2 value of the fit.
 * @result: resulting Chi^2
 * @functionvalue: functionvalue of the fit at a certain point
 * @skips: number of bins that were skipped due to 0 as entry, this needs to be tracked in order to get the ndof
 */
const double FitHandler::getchi2red(){

	double result = 0, functionvalue = 0;
	size_t skips = 0;
	
	//walk over every polynomial
	for(size_t i = 0; i < _m.size(); i++)
	{
		//if the bin entry is 0, it can be a simple misscalculation in Rivet etc.
		if(_b[i] == 0)
		{
			skips++;
			continue;
		}
		
		//calculate the monomials and add them to the result
		for(size_t j = 0; j < _m[0].size(); j++)
			functionvalue += _m[i][j] * _bfp[j];
			
		//get the difference between the functionvalue and the datapoint		
		functionvalue -= _b[i];
		
		//if the uncertainty is 0, the calculation would rise a nan
		//this is prevented by setting it to the arbitrary value of 1e-10
		if(_sigma[i] == 0)
			_sigma[i] = 1e-10;
				
		//calculating the squares of @functionvalue and @sigma, divide them and add this new summand to new overall Chi^2 value @result
		result += (functionvalue * functionvalue) / (_sigma[i] * _sigma[i]);
		
		functionvalue = 0;	
	}
	//return Chi^2 / ndof
	return result / (_m.size() - skips - _m[0].size());
}

/**
 * This function calculates the gradient vectors
 * @pts: container for the anchor points
 * @i: # of the anchor point
 * @functiongradient: noramlized gradient vector
 * @tmp: temporary storage for a monomial the gradient structure
 * @tmpstrc: storage for a new @_structure
 */
vector<double> FitHandler::gradvecfunction(Professor::ParamPoints& pts, size_t i){
	
	//if the static variable @_structure is not initialized yet, its size will be set
	if(pts.structure().size() == 0)
	{
		//the critical pragma prevents multiple accesses to the size at the same time
		#pragma omp critical (structurebuild)
		{
			//set the size to the # of anchor points
			pts.structure().resize(pts.numPoints());
			
			//set for every anchor point the size to the dimension, so there is one for every derivative
			for(size_t j = 0; j < pts.numPoints(); j++)
				pts.structure()[j].resize(pts.dim());
		}
	}
	
	//set up the gradient and resize to the number of dimensions, initialized as 0
	vector<double> functiongradient;
	functiongradient.assign(pts.dim(), 0);
	double tmp = 1;
	
	//calculating the components of the gradient
	for(size_t j = 0; j < functiongradient.size(); j++)
	{
		//ifenough monomials were already calculated, they can be taken directly
		if(pts.structure()[i][j].size() > _bfp.size()) 
			//walking over every monomial and calculate the contribution to the gradient vector
			for(size_t k = 1; k < _bfp.size(); k++)
				functiongradient[j] += _bfp[k] * pts.structure()[i][j][k];
		else
		{
			//if not enough monomials were already calculated in structure, they need to be calculated
			//a critical pragma prevents multiple accesses at the same time and therefore potentially problems
			#pragma omp critical (structureadd)
			{
			//walking over every monomial, calculate what is already available
			for(size_t k = 1; k < pts.structure()[i][j].size(); k++)
				functiongradient[j] += _bfp[k] * pts.structure()[i][j][k];
		
			//use a new variable and calculate the new elements for @_structure there, replace them later
			vector<double> tmpstruc = pts.structure()[i][j];
			
			//set the size to the number of fit parameters as the new maximum needed
			tmpstruc.resize(_bfp.size());
			
			//start the calculation after the last calculated component and walk up to the new maximum needed
			for(size_t k = pts.structure()[i][j].size(); k < tmpstruc.size(); k++)
			{
				//extracting the respective parameter values
				for(size_t l = 0; l < pts.dim(); l++)
					//if the parameter is the one that is derived in the component of the gradient, then its derivative is used
					if(j == l)
						//if the value of the parameter is 0 & the power of the parameter to derive is != 1, the whole monomial will be 0 but the special case in the derivative of 0^0 = 1
						//if the power is 0, the whole monomial will be 0 after derived 
						if((pts.point_scaled(i)[l] == 0 && _power[k][l] != 1) || _power[k][l] == 0)
						{	
							tmp = 0;
							continue;
						}
														
						else			
							//multiply the derived factor
							tmp *= pow(pts.point_scaled(i)[l], _power[k][l] - 1) * _power[k][l];	

					//multiply the rest to @tmp
					else
						//if at least one of the not derived parameters is 0 and its power is !=0, the whole monomial become 0
						if(pts.point_scaled(i)[l] == 0 && _power[k][l] != 0)
						{
							tmp = 0;
							continue;
						}
						//else the rest will be multiplied
						else
							tmp *= pow(pts.point_scaled(i)[l], _power[k][l]);
					
				//put the new monomial part to the @tmpstruc at the specific point in the list
				tmpstruc[k] = tmp;

				//adding the monomial to the gradient component
				functiongradient[j] += tmp * _bfp[k];
				tmp = 1;
			}
			
				//replace the new structure list as the static list for all threads
				//another if-condition concerning the length of the list is meant to prevent possible waiting threads at the beginning of this critical pragma to enter afterwards and replace
				//a longer list by a shorter one
				//therefore this calculation will appear very often in the beginning of the program but fewer and fewer afterwards if more was already calculated. Therefore it is kept static.
				if(tmpstruc.size() > pts.structure()[i][j].size())
					pts.setstructure(i, j, tmpstruc);
			
			}
		}
	}
	
	//normalizing the gradient
	tmp = LinAlg::getabs(functiongradient);
	
	//in order to avoid NaN's, the gradient is only normalized if it isn't a zerovector, else it's normalized vector is the vector itself
	if(tmp != 0)
		for(size_t i = 0; i < functiongradient.size(); i++)
			functiongradient[i] /= tmp;
	
	return functiongradient;
}

/**
 * This function returns Dsmooth
 * @pts: container for the anchor points
 * @gh: container of the gradient vectors of the anchor points
 * @result: result of the calculation that will be returned
 */
const double FitHandler::getDsmooth(Professor::ParamPoints& pts, GradHandler& gh){

	//if the smoothness was not calculated in the iteration ...
	if(_dsmooth_current == 2)
	{
		double result = 0;
		//walking over all anchor points and adding the dotproduct of the gradient vectors
		for(size_t i = 0; i < pts.numPoints(); i++)
			result += LinAlg::dotproduct(gh._gradvecs[i], gradvecfunction(pts, i));
		_dsmooth_current = result / (double) pts.numPoints();
	}
	//returning the mean
	return _dsmooth_current;
}

/**
 * This function calculates the uncertainty of the fit parameters.
 * @num_ipol: number of the bin
 * @mat: storage of the precision matrix, later of the covariance matrix
 * @tmp: temporary storage of precision matrix elements
 */
void FitHandler::setfiterrors(size_t num_ipol){

	//setting up variables
	MatrixXd mat(_bfp.size(), _bfp.size());
	double tmp;

	//fill the matrix; walking over the upper triangle
	for(size_t row = 0; row < _bfp.size(); row++)
		for(size_t col = row; col < _bfp.size(); col++)
		{	
			//pull an element
			tmp = getinvcovmatelement(row, col);
			if(row == col)
				//if the element is on the diagonal then set it
				mat(row, col) = tmp;
			else
			{
				//using the symmetry, setting off-diagonal elements twice
				mat(row, col) = tmp;
				mat(col, row) = tmp;
			}
		}

	//calculating the inverse of the matrix
	mat = mat.fullPivHouseholderQr().inverse();

	//If an element is +/-infinity, it will be set to the maximum of double instead. That way one keeps regular numerical eleements
	for(size_t i = 0; i < _bfp.size(); i++)
		if(mat(i, i) == std::numeric_limits<double>::infinity())
			_bfperr.push_back(std::numeric_limits<double>::max());
		else
			if(-mat(i, i) == std::numeric_limits<double>::infinity())
				_bfperr.push_back(-std::numeric_limits<double>::max());
			else
				_bfperr.push_back(sqrt(mat(i, i)));
				
	//writing the covariance matrix to file
	OutputHandler oh;
	oh.write_covmat(mat, num_ipol);
}

/**
 * This function calculates an element of the precision matrix
 * @i, @j: row/coloumn of the matrix
 * @resulg: numerical value of the matrix element
 */
double FitHandler::getinvcovmatelement(size_t i, size_t j){
	double result = 0;

	//walk over the fit function
	for(size_t k = 0; k < _m.size(); k++)
	{ 
		//If the uncertainty is zero, the result would be inf. Those elements will be skipped. It does not change the overall result, since it will be left out for every element.
		if(_sigma[k] == 0)
			continue;		
		//bin wise calculating d^2\chi^2/da_ida_j with the fit parameters \vec{a}
		result += _m[k][i] * _m[k][j] / (_sigma[k] * _sigma[k]);	
	}
	return result;
}

/**
 * This function calculates the next iterationstep in a bin
 * @pts: container of the anchor points
 * @threshold: threshold for the RR constraint of the fit
 * @kappa: shift in case of applied RR constraint
 */
void FitHandler::nextstep(Professor::ParamPoints& pts, double& threshold, double& kappa){
	
	//reset the smoothness, so that will be calculated again
	_dsmooth_current = 2;
	
	//if every power mentioned in @_power is already in use in @_m, new components of a higher power needs to be calculated
	if(_m[0].size() == _power.size())
	{
		//@_max is the order of the polynom. It will be increased and the new powers will be calculated
		_max++;
		_power = pts.getpower(_max);
	}

	//setting up the new member variables
	increasem(pts);	
	_iterationcounter++;
	makemprime();
	_a.resize(_mprime[0].size());
	expandqr();
	_bprime = LinAlg::normalizevec(_b);
	_bprime = LinAlg::multmatvec(LinAlg::transpose(_q), _bprime);
		
	_d.resize(_m[0].size());
	for(size_t i = 0; i < _d.size(); i++)
		_d[i].assign(_d.size(), 0);
	for(size_t i = 0; i < _d.size(); i++)
		_d[i][i] = 1;
	for(size_t i = 0; i < _a.size(); i++)
	{
		_a[i] = nan("1");
		collinearity(threshold, i, kappa);
	}

	//getting the fit and rescale it
	_bfp = LinAlg::getbestfitparameters(_r, _a, _bprime);	
	rescalebestfitparameters();
}

/**
 * This function walks through an iteration without any checks, so that the matrices are created only
 * @pts: container of the anchor points
 */
void FitHandler::nextstep_walkthrough(Professor::ParamPoints& pts){

	//if every power mentioned in @_power is already in use in @_m, new components of a higher power needs to be calculated
	if(_m[0].size() == _power.size())
	{
		//@_max is the order of the polynom. It will be increased and the new powers will be calculated
		_max++;
		_power = pts.getpower(_max);
	}

	//setting up the new member variables
	increasem(pts);
	_iterationcounter++;
	makemprime();
	expandqr();
}

/**
 * This function sets up the object for a new bin
 * @pts: container of the anchor points
 * @ptvals: bin values for each anchor point
 * @pterrs: error of the bin values for each anchor point
 * @num_ipol: number of the bin
 * @threshold: threshold for the RR constraint of the fit
 * @kappa: shift in case of applied RR constraint
 */
void FitHandler::startbin(Professor::ParamPoints& pts, const vector<double>& ptvals, const vector<double>& pterrs, int& num_ipol, double& threshold, double& kappa){

	//deleting old information
	_iterationcounter = 0;
	_max = 0;
	_bfp.clear();
	_bfperr.clear();
	_m.clear();
	_r.clear();
	_q.clear();
	_d.clear();
	_mprime.clear();
	_a.clear();
	
	//setting up the new reference vector and the corresponding error
	_b = ptvals;
	_sigma = pterrs;
	
	//setting up the power list of zeroth order
	_power = pts.getpower(_max);

	//setting up the matrices and @_a
	increasem(pts);
	makemprime();
	_a.resize(_mprime[0].size());
	expandqr();
	//normalize @_b and preparing it as new righthand side of the LinAlg problem
	_bprime = LinAlg::normalizevec(_b);
	_bprime = LinAlg::multmatvec(LinAlg::transpose(_q), _bprime);

	//setting up the identity matrix
	_d.resize(_m.size());
	for(size_t i = 0; i < _d.size(); i++)
		_d[i].assign(_d.size(), 0);
	for(size_t i = 0; i < _d.size(); i++)
		_d[i][i] = 1;
		
	//checking for RR constraints
	for(size_t i = 0; i < _a.size(); i++)
	{
		_a[i] = nan("1");
		collinearity(threshold, i, kappa);
	}
	//getting the fit and rescaling it
	_bfp = LinAlg::getbestfitparameters(_r, _a, _bprime);
	rescalebestfitparameters();
}

/**
 * This function is a getter for the size of @_bfp
 * @zeros: counter of 0's as fit parameters
 */
size_t FitHandler::getnumfitparams(){
	size_t zeros = 0;
	
	//count the number of 0's at the end of the fit parameters list
	for(size_t i = _bfp.size(); i > 0; i--)
		if(_bfp[i - 1] == 0)
			zeros++;
		else
			break;
		
	//return the number of fit parameters that are != 0	
	return _bfp.size() - zeros;
}

/**
 * This function sets the object up for a certain bin and calculates a given number of iterations
 * @pts: container of the anchor points
 * @ptvals: bin values for each anchor point
 * @pterrs: error of the bin values for each anchor point
 * @num_ipol: number of the bin
 * @threshold: threshold for the RR constraint of the fit
 * @kappa: shift in case of applied RR constraint
 * @bestiteration: number of iterations that will be performed
 */
void FitHandler::setiteration(Professor::ParamPoints& pts, const vector<double>& ptvals, const vector<double>& pterrs, int num_ipol, double threshold, double kappa, int bestiteration){

	//setting up the bin
	startbin(pts, ptvals, pterrs, num_ipol, threshold, kappa);
	
	if(bestiteration > 0)
	{
		//Fast walking the number of iterations along.
		for(size_t i = 0; i < bestiteration - 1; i++)
			nextstep_walkthrough(pts);
		
		//Perform the whole calculation only for the last step. That way the Object is indifferent to a regular stepping.
		nextstep(pts, threshold, kappa);
	}	
}

/**
 * This function calculates the dot product of all normalized gradients at every anchor point with every other
 * @pts: container for the anchor points
 * @gh: container of the gradient vectors of the anchor points
 * @result: vector containing all dot products
 * @tmp: temporary storage of the gradient vectors
 */
vector<double> FitHandler::getallgraddotproducts(Professor::ParamPoints& pts, GradHandler& gh){
	
	vector<double> result;
	vector<vector<double>> tmp;
	tmp.resize(pts.numPoints());
	
	//walk over every possible combination of anchor points
	for(size_t i = 0; i < pts.numPoints(); i++)
	{
		//if a gradient was not calculated yet, do it
		if(tmp[i].empty())
			tmp[i] = gradvecfunction(pts, i);
		for(size_t j = 0; j < pts.numPoints(); j++)
		{
			//if a gradient was not calculated yet, do it
			if(tmp[j].empty())
				tmp[j] = gradvecfunction(pts, j);
			//only calculate the dot product, if the points aren't the same
			if(i != j)
				result.push_back(LinAlg::dotproduct(tmp[i], tmp[j]));
		}
	}
	return result;
}

/**
 * This function adds 0's to the list of fit parameters and their errors in order to fill up a certain order. Afterwards, the parameters will be reordered for Professor 2.2.1 usage.
 * That is needed in order to become usable for Professor 2.2.1.
 * @match_index: mapping list that connects the current order to the order needed by Professor 2.2.1
 * @bfp_sorted, @bfperr_sorted: helper to store the new ordering 
 */
void FitHandler::sortfitparams(vector<int>& match_index){
	
	//push back 0's until there are enough 0's that it matches a certain polynomial order
    for(size_t i = _bfp.size(); i < _power.size(); i++)
		_bfp.push_back(0);
	for(size_t i = _bfperr.size(); i < _power.size(); i++)
		_bfperr.push_back(0);

	//set up the helper
	vector<double> bfp_sorted, bfperr_sorted;
	bfp_sorted.resize(_bfp.size());
	bfperr_sorted.resize(_bfperr.size());
	
	//filling the helper
	for(size_t i = 0; i < match_index.size(); i++)
	{
		bfp_sorted[i] = _bfp[match_index[i]];
		bfperr_sorted[i] = _bfperr[match_index[i]];
	}
	
	//set the ordered list as the entries of the member variables
	_bfp = bfp_sorted;
	_bfperr = bfperr_sorted;
}
