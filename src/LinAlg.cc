#include "Professor/LinAlg.h"

using namespace std;



/**
 * This function returns a coloumn of a matrix
 * @mat: matrix, out of which the coloumn will be extracted
 * @j: number of coloumn
 * @tmp: coloumn that will be returned
 */
vector<double> LinAlg::getcol(vector<vector<double> > mat, size_t j){

	vector<double> tmp;

	for(size_t i = 0; i < mat.size(); i++)
		tmp.push_back(mat[i][j]);

	return tmp;

}

/**
 * This function delivers the absolut value of a vector
 * @vec: vector of interest
 * @abs: absolut value
 */
double LinAlg::getabs(vector<double> vec){

	double abs = 0.;

	//summing up the sqaure of all components of the vector
	for(size_t i = 0; i < vec.size(); i++)
		abs += vec[i] * vec[i];

	//return the sqrt of the sum
	return sqrt(abs);

}

/**
 * This function transposes a given matrix
 * @mat: matrix, that will be transposed
 * @tmp: temporary storage of the transposed matrix
 */
vector<vector<double> > LinAlg::transpose(vector<vector<double> > mat){

	vector<vector<double> > tmp;
	
	//setting @tmp's size to the size of @mat
	tmp.resize(mat[0].size());
	for(size_t i = 0; i < tmp.size(); i++)
		tmp[i].resize(mat.size());

	//transposing by switching the according indices
	for(size_t i = 0; i < mat.size(); i++)
		for(size_t j = 0; j < mat[0].size(); j++)
			tmp[j][i] = mat[i][j];

	return tmp;
}

/**
 * This function multiplicates a matrix and a vector
 * @mat: matrix of the product
 * @vec: vector of the product
 * @result: returned vector, represents the product of @mat and @vec
 * @tmp: helper of summing up a row of @mat and the components of @vec
 */
vector<double> LinAlg::multmatvec(vector<vector<double> > mat, vector<double> vec){

		vector<double> result;
		
		//setting @result's size to the number of rows in @mat
		result.resize(mat.size());
		double tmp = 0;

		//walking over every row of @mat
		for(size_t i = 0; i < mat.size(); i++)
		{
			//Walking over every coloumn of @mat at a given row and walking over every component of @vec. The product is added to @tmp
			for(size_t j = 0; j < vec.size(); j++)
				tmp += mat[i][j] * vec[j];

			//the sum is assigned to one component of the result and tmp is resetted
			result[i] = tmp;
			tmp = 0;
		}

		return result;

}

/**
 * This function normalizes a vector. If the length of the vector is 0, the norm will be the vector itself in order to prevent nan's.
 * @vec: vector, that will be normalized
 * @abs: absolute value of the vector
 * @result: if the length is != 0, the normalized vector will be stored in this variable
 */
vector<double> LinAlg::normalizevec(vector<double> vec){

	//calculate the absolute value of @vec
	double abs = getabs(vec);
	
	//if the length of @vec is != 0, it will be normalized, else the vector itself will be returned 
	if(abs != 0)
	{
		vector<double> result;
		result.resize(vec.size());

		//calculating the normalized components of the vector
		for(size_t i = 0; i < vec.size(); i++)
			result[i] = vec[i] / abs; 
			
		return result;
	}
	else
		return vec;
	

}

/**
 * This function calculates the distance between two vectors
 * @a, @b: vectors of interest
 */
double LinAlg::getdistanceofvectors(vector<double> a, vector<double> b){

	//calculate the difference componentwise
	for(size_t i = 0; i < a.size(); i++)
		a[i] = a[i] - b[i];
		
	//return its absolut value
	return getabs(a);
	
}

/**
 * This function delivers the dotproduct of 2 vectors.
 * The special cases in this function are dotproduct with at least one zerovector. Algebraic, there wouldn't be any problem, but in this program it is.
 * That's because two vectors are exactly the same if they are zerovectors. That is one of the konvergence criteria.
 * Furthermore, if only one of those vectors is a zerovector, this may causes a problem in the context of similarity.
 * @a, @b: vectors for the dotproduct
 * @result: result of the dotproduct
 */
double LinAlg::dotproduct(vector<double> a, vector<double> b){
	
	//if one vector is a zerovector, there will be a notification
	//~ if((getabs(a) == 0 && getabs(b) != 0) || (getabs(a) != 0 && getabs(b) == 0))
		//~ cout << "caution: one zerovector in the dotproduct!" << endl;
	
	//the case of both vectors are zerovectors is handled by explicit definition that the dotproduct of those (normalized) vectors is 1
	if(getabs(a) == 0 && getabs(b) == 0)
		return 1.;
		
	double result = 0;
	
	//summing up the product componentwise
	for(size_t i = 0; i < a.size(); i++)
		result += a[i] * b[i];
	
	return result;

	
}

/**
 * This function solves a problem of the type A*x=b where x is a vector containing the parameters to fit
 * @a: matrix A
 * @x: vector of parameters to fit
 * @b: vector b
 * @result: resulting parameters
 */
vector<double> LinAlg::getbestfitparameters(vector<vector<double> > a, vector<double> x, vector<double> b){

	//setting the size of @result
	vector<double> result;
	result.resize(x.size());

	//starting at the bottom line, walking upwards and calculating the parameters
	for(size_t i = x.size(); i > 0; i--)
	{
		//nan is an indicator in a component of @x for a component, that needs to be calculated
		if(isnan(x[i - 1])){
			//calculating the respective parameter
			result[i - 1] = b[i - 1] / a[i - 1][i - 1];
			
			//forward the solution to every row above the regarded one
			for(size_t j = 0; j < i - 1; j++)
				b[j] -= a[j][i - 1] * result[i - 1];

		}
		else{
			//If a values isn't nan, it was set by @collinearity(). It's value will be forwarded as above.
			result[i - 1] = x[i - 1];
			for(size_t j = 0; j < i - 1; j++)
				b[j] -= a[j][i - 1] * result[i - 1];
		}

	}
	
	return result;

}
