#ifndef __OUTPUTHANDLER__H
#define __OUTPUTHANDLER__H

#include <iostream>
#include <vector>
#include <fstream>
#include <omp.h>
#include "Professor/LinAlg.h"
#include "Professor/ParamPoints.h"
#include "Professor/GradHandler.h"
#include "Professor/FitHandler.h"
#include "Professor/ConfigHandler.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class FitHandler;

/**
* This Class is a container for output functions
*/
class OutputHandler
{
public:

	//Constructor that sets up some necessary parameters
	OutputHandler(const Professor::ParamPoints& pts, int& outdotflag);
	OutputHandler();

	//writes a summary of the fit to the terminal
	void write_binresult(size_t num_ipol, Professor::ParamPoints& pts, GradHandler& gh, FitHandler& fh);

	//writes a dotproduct-summary
	void write_dotproduct(size_t num_ipol, Professor::ParamPoints& pts, GradHandler& gh, FitHandler& fh);

	//writes to the summaryfile
	void write_summary(FitHandler& fh, Professor::ParamPoints& pts, GradHandler& gh);
	
	//setup for the summary file
	void setup_summary();

	//void write_covmat(TMatrixDSym& tmds, size_t num_ipol);
	void write_covmat(MatrixXd& mat, size_t num_ipol);

private:

	/**
	 * @distances: distances of every anchor point to another
	 */
	vector<double> distances;
	
	//calculates the distances between the anchor points
	void setDistances(const Professor::ParamPoints& pts);
};

#endif
