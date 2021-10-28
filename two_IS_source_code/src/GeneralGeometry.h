/*
 * GeneralGeometry.h
 *
 *  Created on: Feb 1, 2017
 *      Author: hornak
 */
//This file delas with simple 3d geometrical tasks



#ifndef GENERALGEOMETRY_H_
#define GENERALGEOMETRY_H_

#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using namespace Eigen;
#include "Cell_parametres.h"
#include <vector>

#include <iostream>
#include <string>
using namespace std;



//Computes the distance between two lines
double distance_two_lines( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector );
//Distance between two segments in 3d
double distance_two_line_segments( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector );
//Same
double distance_two_line_segments( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector , unsigned int interaction_index);
//This function is the same as the previous, but it returns the closest points on the segment and parametres 
double distance_two_line_segments( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector , Vector3d& p_1 , Vector3d& p_2 , double& t_return, double& v_return );
//Distance between the point and a segment
double distance_point_segment( Vector3d a_vector , Vector3d b_vector , Vector3d y_vector );
//Distance between the point and the segment - closest_point_of_segment = closest point on the segment
double distance_point_segment( Vector3d a_vector , Vector3d b_vector , Vector3d y_vector , Vector3d& closest_point_of_segment );
//Distance between the plane and the point
double distance_plane_point( Vector3d plane_axis , Vector3d point_on_plane , Vector3d point );
//Distance between the point and the line(not segment, unlimited line)
double distance_point_line( Vector3d first_line_point , Vector3d second_line_point , Vector3d y_vector );
//Projects the point on the surface 
Vector3d project_point_on_surface_of_elipse( Vector3d position );
//Adds all norms of every bead described by the matrix
double coordinate_norm( MatrixXd coordinates );
//Projects the segment on the plane
Vector3d projection_of_tangent_on_plane( Vector3d plane_axis , Vector3d tangent );
//Projects the point on the plane given by the axis and one point
Vector3d projection_of_point_on_plane( Vector3d plane_axis , Vector3d point_on_plane , Vector3d point_to_be_projected );
//Adds every element in the matrix
double magnitude_3D_matrix( MatrixXd& matrix);



#endif /* GENERALGEOMETRY_H_ */
