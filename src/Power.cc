#include "Professor/Power.h"
#include <iostream>
#include <algorithm>

/**
 * Constructor that sets the constant dimension
 * @size: dimension of the parameters
 */
Power::Power(size_t size): n(size)
{
}

/**
 * This function calculates the list of powers up to a given order
 * and adds them to the overall list @_powerlist
 * @order: order up to which the powers will be calculated
 * @tmp_power: temporary storage for the new list of powers until it will be added to @_powerlist
 * @zero: list of 0's as the 0th order of powers
 * @c: object that delivers the powers for monomials
 */
void Power::setpoweroforder(int order)
{
	
	vector<vector<int>> tmp_power;
	
	//setting the 0th order
	const vector<int> zero(n, 0);
	tmp_power.push_back(zero);
	
	//loop up to the maximum order
	for (int i = 0; i <= order; ++i) 
	{
		//create a @Counter object for the order @i
		Counter c(n, i);
		
		//calculate every combination of powers
		//if they fit the order @i, they will be stored in @tmp_power
		while (c.next(n-1)) 
			if (c.sum() == i) 
				tmp_power.push_back(c.data());     
	}
				
	//saving the calculated list in @_powerlist at the position of the order		
	_powerlist[order] = tmp_power;
}

/**
 * This function is a getter for the list of powers of a given order.
 * If the order isn't calculated yet, it will be calculated within this function.
 * @order: order of the polynom
 */
vector<vector<int>> Power::getpoweroforder(int order)
{
		
	//the length of @_powerlist is equal to the highest order calculated
	//if the requested order is higher, the size of @_powerlist will be adjusted
	if(order > ((int) _powerlist.size()) - 1)
		_powerlist.resize(order + 1);
	
	//if the order is not calculated yet, it will be calculated
	if(_powerlist[order].empty())
		setpoweroforder(order);
	
	return _powerlist[order];
	
}
