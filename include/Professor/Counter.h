#ifndef __COUNTER__H
#define __COUNTER__H

#include <vector>
#include <iostream>
#include <algorithm>


//This code is copied from Professor 2.1.3
//Source: prof2 -> counter.h / counter.cc / Ipol.cc->mkStructure()
using namespace std;

/**
 * This class creates the powers for the parameters in the monomials.
 */
  class Counter {
  public:

    /**
     * This constructor sets up the member variables
     * @dim: this parameter represents the dimension of the parameters
     * @maxval: this is the order of the monomial
     */
    Counter(size_t dim, int maxval) {
      for (unsigned int i=0; i< dim;++i) _data.push_back(0);
      _maxval=maxval;
    };
    
    //Destructor
    ~Counter(void);

	//calculating another step in order to get the powers of a monomial
    bool next(int index);

	//sums up the powers
    int sum();

	//getter for the powers
    vector<int> data() { return _data;}

  private:
	//@_maxval stores the information about the order
	//@_data stores the powers
    int _maxval;
    vector<int> _data;
  };
#endif
