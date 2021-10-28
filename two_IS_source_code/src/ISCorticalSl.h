#ifndef ISCORTICALSL_H_
#define ISCORTICALSL_H_
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseLU>
using Eigen::MatrixXd;
using namespace Eigen;
#include <iostream>
#include <string>
#include "GeneralGeometry.h"
#include "Cell_parametres.h"
#include "KISS.h"
using namespace std;


//there is only one way that microtubule can approach imunological synapse - cortical sliding: from the front
//it corresponds to the fact, that microtubule can be only in the cell

class ISCorticalSl {
public:
    //Default constructor
    ISCorticalSl();
    ISCorticalSl( Vector3d axis_tmp , Vector3d plane_point_tmp , double layer_tmp );
    //This constructor is used
    ISCorticalSl( Vector3d center_front_of_IS_tmp , Vector3d center_rear_of_IS_tmp , double radius_tmp  , double radius_inner_tmp);
    //Copy contructor
    ISCorticalSl(const ISCorticalSl& other);
    
    //Returns the point of the plane
    Vector3d get_point_on_plane();
    //Returns the axis of the IS
    Vector3d get_axis();
    //Returns the center  of the front wall
    Vector3d get_center_front_of_IS();
    //Returns the center  of the rear wall
    Vector3d get_center_rear_of_IS();
    //Gets the inner radius of the IS
    double get_radius_inner();	
    //Returns the radius of the IS
    double get_radius();    
    

    //Overloading perator
    ISCorticalSl& operator=( const ISCorticalSl &other );
    //Computes the distance of the segment in the IS
    double covered_distance_of_segment( Vector3d first_point , Vector3d second_point );
    //Function designed just for visualisation purposes
    std::vector<Vector3d> create_IS_Cortical_Sliding_points( unsigned int number_of_points );  
    //Controls whether it belongs to the IS
    bool check_ISCorticalSl_caught( Vector3d position );    
    

    //Projects the point on the cell membrane
    Vector3d project_point_on_surface( Vector3d position ); 
    //Checks how many of controlled point belongs to the IS
    std::vector<Vector3d> control_IS_Cortical_Sliding_points( std::vector<Vector3d> points_to_control , unsigned int& number_of_IS_points );
    //Checks how many of controlled point belongs to the IS and outside of the IS
    std::vector<Vector3d> control_IS_Cortical_Sliding_points2( std::vector<Vector3d> points_to_control , unsigned int& number_in  , unsigned int& number_out );
    //Simplified control whether the point belongs to the IS
    bool control_IS_Cortical_Sliding_points_for_dynein_testing( Vector3d point );

    //Checks whether point belongs to the cortical sliding IS
    bool control_IS_Cortical_Sliding_points( Vector3d point );
    //Controls the distance from the axis
    bool control_IS_Cortical_Sliding_points_inside_meze( Vector3d point );
    //Function controls whether the points belong to the IS
    std::vector<Vector3d> control_IS_Cortical_Sliding_points_2_b( std::vector<Vector3d> points_to_control , std::vector<Vector3d>& inside_points );
    //Controls whether points belong to the IS
    std::vector<Vector3d> control_IS_Cortical_Sliding_points_2_b_for_dynein_testing( std::vector<Vector3d> points_to_control , std::vector<Vector3d>& inside_points );    
    //Controls whether the points are in the IS or at the margins
    std::vector<Vector3d> control_of_inside_and_margins( std::vector<Vector3d> points_to_control , std::vector<Vector3d>& points_margins );
    

	virtual ~ISCorticalSl();
private:
    Vector3d axis;
    Vector3d plane_point;
    double layer;
    double radius;
    double radius_inner;
    Vector3d center_front_of_IS;
    Vector3d center_rear_of_IS;	


};

#endif /* ISCORTICALSL_H_ */
