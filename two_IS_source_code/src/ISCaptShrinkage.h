/*
 * ISPlanar.h
 *
 *  Created on: Oct 10, 2016
 *      Author: hornak
 */

#ifndef ISCAPTSHRINKAGE_H_
#define ISCAPTSHRINKAGE_H_
#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using namespace Eigen;
#include <iostream>
#include <string>
#include "GeneralGeometry.h"

using namespace std;



class IS_Capt_Shrinkage {
public:
        //Default constructor
	IS_Capt_Shrinkage();
	//The cylinder of the capture-shrinkage IS is defined by the centers of the front, rear side and the radius 
	IS_Capt_Shrinkage( Vector3d center_of_IS_front , double radius_argument , Vector3d center_of_IS_rear );
	IS_Capt_Shrinkage( Vector3d center_of_IS_front , double radius_argument , Vector3d center_of_IS_rear , Vector3d trap_arg );
	//Copy constructor 
	IS_Capt_Shrinkage( IS_Capt_Shrinkage& tmp );
	//Overloading the operator
	IS_Capt_Shrinkage operator=( const IS_Capt_Shrinkage& tmp );
	//Gets the center of the front side of the IS
	Vector3d get_center_of_IS_front() const;
	//Gets the center of the rear side of the IS
	Vector3d get_center_of_IS_rear() const;
	//Gets the radius of the IS
	double get_radius_of_IS() const;
	//Gets the axis of the IS
	Vector3d get_axis_of_IS() const;
	//Sets the center of the front side
	void set_center_of_IS_front( Vector3d center_argument_front );
	//Sets the center of the rear side	
	void set_center_of_IS_rear( Vector3d center_argument_rear );
        //Sets the radius of the IS
	void set_radius_of_IS( double radius_center );
	//Sets the axis of the IS
	void set_axis_of_IS( Vector3d axis_argument );
	//Checks whether the segment intersect the IS
	bool check_IS_capture_shrinkage_caught_Main( Vector3d position );
	//Check whether the segment lies in the proximity of the IS
        bool check_IS_capture_segment_control( Vector3d first_point , Vector3d tangent );
        //Destructor
	virtual ~IS_Capt_Shrinkage();

private:

	//The IS capture-shrinkage has a shape of a cylinder 
	//Dyneins are located at the intersection of the cylinder and the sphere of plasma membrane 
	//The cylinder is given by the centers of the front(close to the IS) and the rear side
	//by the axis and the radius
	Vector3d center_of_IS_front;
	Vector3d center_of_IS_rear;
	Vector3d axis_of_IS;
	double radius_of_IS;
	


};

#endif 
