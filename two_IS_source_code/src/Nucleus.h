/*
 * Nucleus.h
 *
 *  Created on: Jan 18, 2017
 *      Author: hornak
 */

#ifndef NUCLEUS_H_
#define NUCLEUS_H_
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include "Nucleus_parametres.h"
#include "Cell_parametres.h"


using Eigen::MatrixXd;
using namespace Eigen;

class Nucleus {
public:
	//Default contructor
	Nucleus();
	//Constructor used in the simulations
	Nucleus( Vector3d center_argument );
	//Copy contructor
	Nucleus& operator=( const Nucleus &tmp );
	//Gets the axis of the nucleus
	double get_A_Axis() const;
	//Gets the axis of the nucleus
	double get_B_Axis() const;
	//Gets the center
	Vector3d get_center();
	//The force acting on a point of the bead of the microtubule
	Vector3d force_position_nucleus_wall( Vector3d position );
	//The force acting on a point of the bead of the MTOC
        Vector3d force_position_nucleus_wall_MTOC( Vector3d position );


	virtual ~Nucleus();

private:
	//Nucleus is an ellipsoid defined by two axis and the center
	double a_axis;
	double b_axis;
	Vector3d center;
};

#endif /* NUCLEUS_H_ */
