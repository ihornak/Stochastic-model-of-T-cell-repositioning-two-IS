/*
 * CellShape.h
 *
 *  Created on: Apr 28, 2017
 *      Author: hornak
 */

#ifndef CELLSHAPE_H_
#define CELLSHAPE_H_
#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h>

using namespace std;


//Cell shape plays a role only if cell is not a sphere
class Cell_Shape {
public:
	Cell_Shape();
	Cell_Shape( double A_Axis_arg , double B_Axis_arg , double cylinder_width_arg );

	double calculate_z_cylinder();
	double get_A_Axis();
	double get_B_Axis();
	double get_cylinder_width();
	double get_z_value();

	virtual ~Cell_Shape();
private:
	double A_Axis;
	double B_Axis;
	double cylinder_width;
	double z_value;

};

#endif /* CELLSHAPE_H_ */
