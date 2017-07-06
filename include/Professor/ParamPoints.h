#ifndef PROF_PARAMPOINTS_H
#define PROF_PARAMPOINTS_H

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>
#include <stdexcept>
#include "Professor/LinAlg.h"
#include "Professor/Power.h"

namespace Professor {


  /// Throwable error
  struct ParamPointsError : public std::runtime_error {
    ParamPointsError(const std::string& reason) : std::runtime_error(reason) { }
  };


  /// Typedef for a list of parameters, defining a parameter point
  typedef std::vector<double> ParamPoint;


  /// @todo I think we need to rename the more structured object
  // typedef const std::vector< std::vector<double> > ParamPointVec;


  /// Class for the parametrisation hypercube, i.e. anchors
  class ParamPoints {
  public:

    /// @todo Also record the parameter names in this object

    /// Constructor, takes anchors as vector<vector<double>> = vector<ParamPoint>
    /// @todo Also take the parameter names as an arg
    ParamPoints(const std::vector< std::vector<double> >& ppoints);

    /// Implicit conversion operator to vector<vector<double>> = vector<ParamPoint>
    operator std::vector< std::vector<double> > () const {
      return _parampoints;
    }


    /// Overly complicated push_back
    void addParamPoint(const std::vector<double>& p) {
      if (_locked) {
        /// @todo Throw a ParamPointsError or similar rather than calling abort()
        std::cerr << "Adding point to locked collection not implemented, aborting" << std::endl;
        abort();
      }
      // This ensures that all ppoints are of the same dimension
      if (_parampoints.size() > 0)
        assert(p.size() == _parampoints[0].size()); ///< @todo Throw an exception instead of assert
      _parampoints.push_back(p);
    }

    /// Number of anchor points
    int numPoints() const { return _parampoints.size(); }

    /// Dimension of (anchor) points
    int dim() const {
      assert(!_parampoints.empty()); //< Emptiness should not be possible
      return _parampoints.front().size();
    }

    /// Centre of the anchor hyper cube
    std::vector<double> ptcenters() const;

    /// Lowest edge of anchor hyper cube
    std::vector<double> ptmins() const;

    /// Top edge of anchor hyper cube
    std::vector<double> ptmaxs() const;

    /// Edges of the anchor hyper cube
    std::vector< std::pair<double, double> > ptedges() const;

    /// print message: anchors
    /// @todo These non-redirectable print functions are a bad idea
    void printPoints() const;

    /// print message: meta info
    /// @todo These non-redirectable print functions are a bad idea
    void printMeta() const;

    /// @todo Add toString method(s) for the main logic in the print functions above, and to connect to Python's ParamPoints.__str__ function.

    /// Remove all anchors
    void reset() { _parampoints.clear(); _locked = false; };

    /// Header representation with metadata necesary to write out ProfDF
    std::string toString(const std::string& info="") const {
      std::stringstream ss;
      if (!info.empty()) ss << "# INFO " << info << "\n";
      ss << "# MINV ";
      for (const double& a : ptmins()) ss << a<< " ";
      ss << " \n";
      ss << "# MAXV ";
      for (const double& a : ptmaxs()) ss << a<< " ";
      ss << " \n";
      // Slightly redundant but nice to have as consistency check
      ss << "# DIM ";
      ss << dim();
      ss << "\n";
      // Write out parameter names if defined
      if (names().size()>0) {
        ss << "# PARAMS ";
        for (const std::string& s : names()) ss << s<< " ";
        ss << "\n";
      }

      return ss.str();
    }

    /// Get all anchor points
    const std::vector< std::vector<double> >& points() const { return _parampoints; }

    /// Get parameter names if set
    const std::vector< std::string >& names() const { return _names; }

    void setNames(std::vector<std::string >);

    /// Get one anchor point
    /// @todo Generalise to other sorts of key lookup? Needed? Bounds checking / checking key existence
    const std::vector<double>& point(size_t i) const { return points().at(i); }

    // Setter for the indices referring to a hypercube surrounding a point with index @i
    std::vector< size_t>& gethypercube(size_t i);

    // Getter of @_hypercubes
    std::vector< std::vector<size_t>>& hypercubes() {return _hypercubes;}

    // Getter of @_structure
    std::vector< std::vector< std::vector<double>>>& structure() {return _structure;}

    // This function adds vectors to @_structure
    void setstructure(size_t i, size_t j, std::vector<double> vec);

    // Getter of @_parampoints_scaled
    const std::vector<std::vector<double>>& points_scaled() const { return _parampoints_scaled;}

    // Getter of point @i out of @_parampoints_scaled
    const std::vector<double>& point_scaled(size_t i) const { return points_scaled().at(i); }
    
    // This function rescales @_parampoints and stores them in @_parampoints_scaled
    void rescale();

    // This function gets the list of powers for every parameter in every term of the polynomial function at given order @order
    std::vector<std::vector<int>> getpower(int order) {return _pow.getpoweroforder(order);}

    // Getter of @_pow
    Power& getpow() {return _pow;}

    void clearall();

  private:

    // structure of the powers of the gradient vectors of the polynomial function
    std::vector< std::vector< std::vector<double>>> _structure;

    // list of indices for hypercubes of different points
    std::vector< std::vector<size_t>> _hypercubes;

    // parameter and scaled parameter values
    std::vector< std::vector<double> > _parampoints, _parampoints_scaled;

    std::vector< std::string > _names;

    bool _locked;

    // Helper class to retrieve the list of powers
    Power _pow;

  };


}

#endif
