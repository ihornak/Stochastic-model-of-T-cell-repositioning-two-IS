/*
 * MTOC2.h
 *
 *  Created on: May 30, 2017
 *      Author: hornak
 */

#ifndef MTOC2_H_
#define MTOC2_H_


#include <iostream>
#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using namespace Eigen;
#include <vector>
#include <random>
#include <iostream>



#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseLU>
using namespace std;

#include "simulationofCell.h"
#include "MTOCparam.h"
#include "GeneralGeometry.h"
#include "mersenne.h"
#include "simulationofCell.h"
//#include <eigen3/Eigen/src/SparseLU>





class MTOC2
{
public:
    //Default constructor
    MTOC2();
    MTOC2( unsigned int number_of_points_arg , Vector3d center_of_MTOC , string plane );

    //Overloading the operator	
    MTOC2& operator=( const MTOC2& tmp );
	
    //Setters and getters
    //Sets the coordinates of the MTOC
    void set_coordinates( MatrixXd& coordinates_arg );
    //Returns the number of sprouting points
    unsigned int get_number_of_points() const;
    //Returns the coordinates of all MTOC points
    MatrixXd get_coordinates() const;
    //Returns number-th point of the MTOC
    Vector3d get_point( unsigned int number ) const;
    //Returns the center of the MTOC
    Vector3d get_center() const;
    //Returns the radius of the MTOC
    double get_radius() const;
    //Sets the radius of the MTOC    
    double set_radius(  );   
    //Sets the time
    void set_time_clock( double time_tmp ){ this->time_clock = time_tmp;  };     
    //Gets the point of the original coordinates of the MTOC
    Vector3d get_point_original( unsigned int number ) const;
    //Gets the original orientation of the sprouting point from the MTOC
    Vector3d get_original_orientation( unsigned int number ) const;
    //Returns the original center of the MTOC
    Vector3d get_center_original( ) const;
    //Returns the original coordinates of the MTOC at the beginning of the simulation
    MatrixXd get_original_coordinates() const;
    //This returns the average position of the MTOC beads - used only for the control
    Vector3d get_center_of_gravity() const;
    //Returns the effective friction of the MTOC beads 
    double get_effectiveFriction() const { return this->effectiveFriction; }
    //Sets the friction of one microtubule bead
    void set_friction();
    //This will be used for 3d MTOC
    Vector3d get_side_center( unsigned int side );
    //Gets the exactly opposite point of the MTOC
    unsigned int get_opposite_point_index( unsigned int point_index );
    //Gets the point at the opposite side of the MTOC - opposite but parallel to the axis
    unsigned int get_opposite_direct_point_index( unsigned int point_index );

    //Returns the tangent connecting the center and the point
    Vector3d tangent_to_center( unsigned int number );
    //Gets the tangent connecting the bead to the one with higher index
    Vector3d tangent_from_next_bead( unsigned int number );
    //Returns the tangent between the two beads
    Vector3d tangent_between_beads(  unsigned int number_1 , unsigned int number_2 );
    //Controls the extent of numerical impressisions 
    void controlMTOC( );
    //Basic resizing of tangent lengths
    void adjustMTOC( );
    //Gets the random point from the opposite side of the MTOC
    unsigned int get_random_point_on_opposite_side( unsigned int point_index );
    //Returns a random point chosen from the five most distant points of the MTOC
    unsigned int get_random_point_on_opposite_side_without_bias( unsigned int point_index );
    //Comptation of the projection matrixes 
    void getMatrixes_Sparse( MatrixXd &projection_Matrix );
    void getMatrixes_Sparse_2( MatrixXd &projection_Matrix );
    void getMatrixes_Sparse_3( MatrixXd &projection_Matrix );
    void getMatrixes_Sparse_4( MatrixXd &projection_Matrix );
    //Get bending forces acting on the MTOC beads
    void getRandomForces(  double timeStep , MatrixXd &randomForces );
    //Gets the average distance from the plane parallel to the plane of the MTOC intersecting the center of coordinates - for control
    double get_average_altitude();
    //Sum of norms of all forces acting on the beads - just for control
    double magnitude_of_the_force( MatrixXd& Forces  );

    MatrixXd rotate_MTOC_from_original_orientations();
    MatrixXd rotate_MTOC_from_original_orientations_2();


    //Prints all tangents to the terminal - basic control
    void print_norm_tangent_to_center();
    //Prints all points to the terminal - basic control
    void print_points();


    //One ste using Euler algorithm - just for testing
    void Euler_algorithm( MatrixXd& Forces );
    //First half of the Mid-Step algorithm
    void oneStepMidStepAlgorithm_1_half( MatrixXd& Forces );
    //Second half of the Mid-Step algorithm
    void oneStepMidStepAlgorithm_2_half( MatrixXd& Forces , MatrixXd& coordinates_arg);
    //One step using Mid-Step algorithm
    void oneStepMidStepAlgorithm( );
    //Propagation of the isolated MTOC - just for the control
    void MidStepAlgorithm( double Time_0 , double Time_1 );



    //Resize the MTOC due to the numerical imprecisions 
    void resize_from_originals();
    //Rezizes the MTOC and keeps the center on the axis
    void resize_from_originals_axis_calibration();






    virtual ~MTOC2();
	
//Currently unused costructors
    MTOC2( unsigned int number_of_points_arg );
    MTOC2( unsigned int number_of_points_arg , Vector3d center_of_MTOC  );    
    MTOC2( unsigned int number_of_points_arg , Vector3d center_of_MTOC , string plane  , string plane_2 );	
    //So far, not used  - it is for later for 3d MTOC	
    double set_polygon_angle(  );
    //So far, not used  - it is for later for 3d MTOC	    
    void set_axis_of_rotation();	
	

private:
	unsigned int number_of_points;
	MatrixXd coordinates;
	double radius_of_MTOC;
  	double effectiveFriction;
  	//original coordinates of the MTOC
  	MatrixXd original_orientations;
  	//Original orientation of the MTOC
  	Vector3d original_MTOC_orientation;
  	std::vector<Vector3d> axis_of_rotation;
	MatrixXd original_coordinates;
	//So far, not used  - it is for later for 3d MTOC
	double polygon_angle = 3.1415926535897 / 2.0 / 4.0;
	//Time of the simulation	
	double time_clock;


};

#endif /* MTOC2_H_ */
