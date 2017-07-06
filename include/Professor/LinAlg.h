#ifndef __LINALG__H
#define __LINALG__H

#include <iostream>
//#include <stdlib.h>
#include <vector>
#include <math.h>


using namespace std;

/**
* This Class is a container for linear algebra.
*/
class LinAlg
{
public:

	//extracts a coloumn of a matrix
	static vector<double> getcol(vector<vector<double> > mat, size_t j);

	//calculates the absolut value of a vector
	static double getabs(vector<double> vec);

	//transposes a matrix
	static vector<vector<double> > transpose(vector<vector<double> > mat);

	//multiplicates a matrix and a vector
	static vector<double> multmatvec(vector<vector<double> > mat, vector<double> vec);

	//normalizes a vector
	static vector<double> normalizevec(vector<double> vec);

	//calculates the difference between two vectors
	static double getdistanceofvectors(vector<double> a, vector<double> b);

	//calculates the dotproduct of two vectors
	static double dotproduct(vector<double> a, vector<double> b);

	//solves a problem of the type matrix * vector = vector
	static vector<double> getbestfitparameters(vector<vector<double> > a, vector<double> x, vector<double> b);

};

#endif
