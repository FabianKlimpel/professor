#include "Professor/ParamPoints.h"
#include <algorithm>

namespace Professor {

  using namespace std;


  ParamPoints::ParamPoints(const vector< vector<double> >& ppoints) {
    /// @todo Throw a ParamPointsError or similar rather than this assert
    assert(!ppoints.empty());
    // _parampoints.clear();
    // _locked = false;
    // for (size_t i = 0; i < p.size(); ++i) {
    //   _parampoints.push_back(p[i]);
    // }
    _parampoints = ppoints;
    _locked = true;
    _pow.setdim(dim());
  }

  vector<double> ParamPoints::ptcenters() const {
    vector<double> temp_max, temp_min;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      vector<double> temp;
      for (size_t j = 0; j < numPoints(); j++) { // iteration over anchors
        temp.push_back(_parampoints[j][i]);
      }
      temp_max.push_back(*max_element(temp.begin(), temp.end()));
      temp_min.push_back(*min_element(temp.begin(), temp.end()));
    }
    vector<double> center;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      center.push_back(temp_min[i] + 0.5* (temp_max[i] - temp_min[i]));
    }

    return center;
  }


  vector<double> ParamPoints::ptmins() const {
    vector<double> temp_min;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      vector<double> temp;
      for (size_t j = 0; j < numPoints(); j++) { // iteration over anchors
        temp.push_back(_parampoints[j][i]);
      }
      temp_min.push_back(*min_element(temp.begin(), temp.end()));
    }
    return temp_min;
  }


  vector<double> ParamPoints::ptmaxs() const {
    vector<double> temp_max;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      vector<double> temp;
      for (size_t j = 0; j < numPoints(); j++) { // iteration over anchors
        temp.push_back(_parampoints[j][i]);
      }
      temp_max.push_back(*max_element(temp.begin(), temp.end()));
    }
    return temp_max;
  }


  vector< pair<double, double> > ParamPoints::ptedges() const {
    const vector<double> mins = ptmins();
    const vector<double> maxs = ptmaxs();
    vector< pair<double, double> > edge;
    edge.reserve(dim());
    for (size_t i = 0; i < dim(); i++) {
      edge.push_back( pair<double, double>(mins[i], maxs[i]));
    }
    return edge;
  }

  void ParamPoints::setNames(std::vector<std::string > names) {
    if (_names.size() == 0) { // No names set so far
      if (dim() == names.size()) { // Sanity check
        for (size_t i = 0; i < names.size(); i++) {
          _names.push_back(names[i]);
        }
      }
      else {
        stringstream ss;
        ss << "ParamPoints::setNames: dimension mismatch (" << dim() << "dimensions vs. " << names.size() << " names)  ";
        throw ParamPointsError(ss.str());
      }
    }
    else { // Names already set
      stringstream ss;
      ss << "ParamPoints::setNames: Names already set!";
      throw ParamPointsError(ss.str());
    }
  }

  void ParamPoints::printMeta() const {
    cout << "Nr. of points: " << numPoints() << endl;
    cout << "Dimension:     " << dim() << endl;
    cout << "Center:       ";
    for (size_t i = 0; i < dim(); i++) {
      cout << " " << ptcenters()[i];
    }
    cout << endl;
    cout << "Edges:" << endl;
    const vector< pair<double, double> > edges = ptedges();
    for (size_t i = 0; i < dim(); i++) {
      cout << edges[i].first << " < " << edges[i].second << endl;
    }
    cout << endl;
  }


  void ParamPoints::printPoints() const {
    for (size_t i = 0; i < numPoints(); ++i) {
      cout << "Point " << i << ":" << endl;
      for (size_t j = 0; j < dim(); ++j) {
        cout << _parampoints[i][j] << " ";
      }
      cout << endl;
    }
  }
  
  /**
   * This function calculates a list of indices referring to the hypercube surrounding the point with index @i
   * @i: number of the point around which a hypercube should be constructed
   * @center: temporary storage of the coordinates of point @i
   * @result: storage of the hypercube indices
   * @distances: list of the shortest distances found in every direction
   * @bin: this parameter is used as an indicator for the direction of the proposal point in a binary coding
   */
  vector<size_t>& ParamPoints::gethypercube(size_t i){
	
	//If the size of the hypercubes does not fit, it will be resized.
	//This construction is meant as initializer of @_hypercubes
	if(_hypercubes.size() != _parampoints.size())
		_hypercubes.resize(_parampoints.size());
		
	//if the hypercube is already calculated, return it
	if(!_hypercubes[i].empty())
		return _hypercubes[i];
	
	//setting up @center and the resultvector
	vector<double> center = _parampoints[i];
	vector<size_t> result;
	
	//setting the resultsize as 2 * #dimensions, so that for every dimension at least 2 points cann be selected
	result.resize(pow(2., _parampoints[0].size()));
	
	//if a point isn't set, it's value is #number of dimensions + 1, so that these points can be found and won't be used for further calculations
	//this is important for points at the border of the sampleregion
	result.assign(result.size(), _parampoints.size());
	
	//The same procedure as before is done for the distances, so that a check is possible if the smallest possible hypercube cann be constructed. Therefore it's initial values are set to infinity.
	vector<double> distances;
	distances.resize(pow(2., _parampoints[0].size()));
	distances.assign(distances.size(), numeric_limits<double>::infinity());	
	size_t bin = 0;
	
	for(size_t j = 0; j < _parampoints.size(); j++)
		//if a datapoint matches the @center, it will be skipped
		if(i == j)
			continue;
		else
		{
			//component wise check, if the components are bigger than the center point
			for(size_t k = 0; k < _parampoints[0].size(); k++)
				if(center[k] > _parampoints[j][k])
					//All 'if' checks were binary questions, therefore a binary representation as summary can be found.
					//@bin is a value that can directly connect to the overall situation of the test point due to its value
					bin += pow(2., k);
		
			//check up, if the new trial point is closer than the stored one in the category specified by @bin
			if(LinAlg::getdistanceofvectors(center, _parampoints[j]) < distances[bin])
			{
				//if closer than set it as new point
				distances[bin] = LinAlg::getdistanceofvectors(center, _parampoints[j]);
				result[bin] = j;				
			}
			
			bin = 0;
		}
		
	//write in a critical enivorment, so that it is thread safe
	#pragma omp critical (sethypercube)
		//store the result
		_hypercubes[i] = result;
	
	return _hypercubes[i];
	
}
	
/**
 * This function adds a new vector to the list of structures
 * @i, @j: indices to locate the storage position
 * @vec: vector that will be added
 */
void ParamPoints::setstructure(size_t i, size_t j, vector<double> vec){_structure[i][j] = vec;}
	
/**
 * This function maps the parameter points @_parampoints onto a [0,1]-hypercube and store them in @_parampoints_scaled
 * @diff: temporary storage of the difference between the maximum and minimum in every dimension
 */
void ParamPoints::rescale(){
	
	//if the points weren't calculated yet, do it
	if(points_scaled().empty())
	{
		//resize @_parampoints_scaled so that it can be easy set
		_parampoints_scaled.resize(numPoints());
		for(size_t i = 0; i < numPoints(); i++)
			_parampoints_scaled[i].resize(dim());
		
		//walk over every dimension and calculate the difference between the maximum and the minimum value
		vector<double> diff;
		for(size_t i = 0; i < dim(); i++)
			diff.push_back(ptmaxs()[i] - ptmins()[i]);
		
		//walk over every parameter set dimensionwise and rescale it onto a [0,1]-hypercube
		for(size_t i = 0; i < numPoints(); i++)
			for(size_t j = 0; j < dim(); j++)
				_parampoints_scaled[i][j] = (points()[i][j] - ptmins()[j]) / diff[j];	
	}		
}

void ParamPoints::clearall(){

	_structure.clear();
	_hypercubes.clear();
	_parampoints.clear();
	_parampoints_scaled.clear();
	_names.clear();
	_pow.clearall();
}
}
