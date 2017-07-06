#include "Professor/OutputHandler.h"

using namespace std;

/**
 * Constructor
 * Here, the sizes of the vectors that will later contain fitparameters etc. will be set. Additionally, the distances and centers will be calculated
 * @ch: carries the flag for the dotproduct-summary
 * @rh: handler of the reference data
 * @minparval: minimum of the parametervalues
 * @maxparval: maximum of the parametervalues
 */
OutputHandler::OutputHandler(const Professor::ParamPoints& pts, int& outdotflag){	
	//if the dotproduct-summary should be written, the distances will be calculated
	if(outdotflag)
	{
		//Resizing the @distances vectors. No distances between an anchor point and itself will be calculated, therefore the size is adapted accordingly.
		distances.resize(pts.numPoints() * pts.numPoints() - pts.numPoints());
		//Calculate the distances
		setDistances(pts);
	}
}	

OutputHandler::OutputHandler(){}

/**
 * This function calculates the distances between the anchors points
 * @ch: carries the number of MC runs
 * @rh: handles the reference data
 * @skips: no distances between an anchor point and itself is calculated, so the index shift needs to be registered by this variable
 */
void OutputHandler::setDistances(const Professor::ParamPoints& pts){
	//walk over every point in parallel
	#pragma omp parallel for schedule(dynamic)
	for(size_t i = 0; i < pts.numPoints(); i++)
	{		
		size_t skips = 0;
		//compare every point with every other
		for(size_t j = 0; j < pts.numPoints(); j++)
			//if both points are the same, the iteration is skipped
			if(i == j)
				skips--;				
			else
				//calculate the distance between the vectors and store them
				distances[i * pts.numPoints() - i - 1 + j + skips] = LinAlg::getdistanceofvectors(pts.point_scaled(i), pts.point_scaled(j));
	}
}

/**
 * This function creates the summary file and writes its header
 * @outsummary: delivers the output to file functionality
 */
void OutputHandler::setup_summary(){	
	
	string line;
	ifstream insummary;
	insummary.open("summary");
	getline(insummary, line);
	insummary.close();
	if(line.substr(0, 3) == "Chi")
		return;

	//set up the file and write the header
	ofstream outsummary;
	outsummary.open("summary");
	outsummary << "Chi^2" << "\t" << "Chi2^2,red" << "\t" << "Iterations" << "\t" << "Dsmooth" << endl;
	outsummary.close();
}

/**
 * This function write the result of the fit for a bin to the console
 * @num_ipol: # of the bin
 * @ch: container for the observable name and the analysis name
 * @rh: handler of the reference data
 * @fh: handler of the fit
 * @store: counter to identify the right analysis/observable
 */
void OutputHandler::write_binresult(size_t num_ipol, Professor::ParamPoints& pts, GradHandler& gh, FitHandler& fh){
	cout << endl << "Result for bin " << num_ipol << ":" << endl;

	//the constrained monomial numbers will be printed to the terminal
	cout << "RR constraint:\t";
	for(size_t i = 0; i < fh.getnumfitparams(); i++)
		if(!std::isnan(fh.geta()[i]))
			cout << i << "\t";
	cout << endl;
	
	//write further summary variables will be printed to the terminal
	cout << "Dsmooth:\t\t" << fh.getDsmooth(pts, gh) << endl;
	cout << "chi2:\t\t\t" << fh.getchi2() << endl;
	cout << "iterationcounter:\t" << fh.getiterationcounter() << endl;
	cout << "max. power:\t\t" << fh.getmax() << endl;
	cout << "-------------------------" << endl;	
}

/**
 * This function writes the dotproduct-summary
 * @num_ipol: # of the bin
 * @ch: contains the number of runs
 * @rh: handles the reference data
 * @fh: handles the fit
 * @outdot: delivers the output to file functionality
 * @rhdot, @fhdot: stores the dot products of the reference data and the fit at every anchor point
 */
void OutputHandler::write_dotproduct(size_t num_ipol, Professor::ParamPoints& pts, GradHandler& gh, FitHandler& fh){
	//set up the output
	ofstream outdot;
	outdot.open(("dotproduct" + to_string(num_ipol)).c_str());
	
	//calculate all dot products
	vector<double> rhdot = gh.getallgraddotproducts(), fhdot = fh.getallgraddotproducts(pts, gh);
	
	//write some summary parameters
	outdot << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << fh.getDsmooth(pts, gh) << "\t" << fh.getchi2() << "\t" 
		<< fh.getchi2red() << "\t" << fh.getiterationcounter() << endl;
		
	//write the dot products of the reference data, the fit and the distances between the anchor points used for the dot product
	for(size_t k = 0; k < rhdot.size(); k++)
		outdot << rhdot[k] << "\t" << fhdot[k] << "\t" << distances[k] << "\t";
	outdot << endl;
						
	outdot.close();
}
	
/**
 * This function writes a summary of a fit to the summary file
 * @fh: handles the fit
 * @rh: handles the reference data
 * @outsummary: delivers the output to file functionality
 */
void OutputHandler::write_summary(FitHandler& fh, Professor::ParamPoints& pts, GradHandler& gh){
	//open the file and continue writing at its end
	ofstream outsummary;
	outsummary.open("summary", ofstream::out | ofstream::app);
	//write out the resulting Chi^2, the reduced Chi^2, the number of iterations needed and the smoothness
	outsummary << fh.getchi2() << "\t" << fh.getchi2red() << "\t" << fh.getiterationcounter() << "\t" << fh.getDsmooth(pts, gh) << endl;
	outsummary.close();
}

void OutputHandler::write_covmat(MatrixXd& mat, size_t num_ipol){
	ofstream outcovmat;
	outcovmat.open(("covmat_" + to_string(num_ipol)).c_str());
	
	for(size_t row = 0; row < (size_t) mat.rows(); row++)
	{
		for(size_t col = 0; col < (size_t) mat.cols(); col++)
			if(mat(row, col) == std::numeric_limits<double>::infinity())
				outcovmat << std::numeric_limits<double>::max() << " ";
			else
				if(-mat(row, col) == std::numeric_limits<double>::infinity())
					outcovmat << -std::numeric_limits<double>::max() << " ";
				else
					outcovmat << mat(row, col) << " ";
		outcovmat << "\n";
	}
	
	outcovmat.close();
}

//~ void OutputHandler::write_covmat(TMatrixDSym& tmds, size_t num_ipol){
	//~ ofstream outcovmat;
	//~ outcovmat.open(("covmat_" + to_string(num_ipol)).c_str());
	//~ 
	//~ for(size_t row = 0; row < (size_t) tmds.GetNrows(); row++)
	//~ {
		//~ for(size_t col = 0; col < (size_t) tmds.GetNcols(); col++)
			//~ if(tmds[row][col] == std::numeric_limits<double>::infinity())
				//~ outcovmat << std::numeric_limits<double>::max() << " ";
			//~ else
				//~ if(-tmds[row][col] == std::numeric_limits<double>::infinity())
					//~ outcovmat << -std::numeric_limits<double>::max() << " ";
				//~ else
					//~ outcovmat << tmds[row][col] << " ";
		//~ outcovmat << "\n";
	//~ }
	//~ 
	//~ outcovmat.close();
	//~ 
//~ }
