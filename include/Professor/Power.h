#ifndef __POWER__H
#define __POWER__H


#include <vector>
#include "Professor/Counter.h"
#include <iostream>

using namespace std;

/**
 * This class serves as a look up table for the list of powers of a certain order
 */
  class Power {
  public:

   	//Constructor
	Power(size_t size);
	Power() {};

	//getter for the list of powers of a certain order
	vector<vector<int>> getpoweroforder(int order);

	vector<vector<vector<int>>>& getpowerlist() {return _powerlist;}

   	void setdim(size_t size) {n = size;}

	void clearall() {_powerlist.clear();}
  private:

	//setter for the list of powers of a new order
	void setpoweroforder(int order);

	//list of lists of powers of certain orders
	vector<vector<vector<int>>> _powerlist;
	
	//dimension of the polynom
	size_t n;

  };
#endif
