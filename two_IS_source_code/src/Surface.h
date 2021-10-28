/*
 * Copyright 2017 <copyright holder> <email>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SURFACE_H
#define SURFACE_H


#include "simulationofCell.h"
#include "Cell_parametres.h"
#include "KISS.h"
#include "IS_Dynein_Cell_surface.h"
#include "IS_Cortical_Sl_parameter.h"
#include "GeneralGeometry.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include "ISCaptShrinkage.h"
#include "ISCorticalSl.h"



using Eigen::MatrixXd;
using namespace Eigen;
#include "Dynein.h"

using namespace std;


class Surface
{
public:
//Default constuctor
Surface();
//This is the constructor for capture-shrinkage IS
Surface(  double density , double width , IS_Capt_Shrinkage IS_tmp );
//Constructor for cortical sliding dyneins on the surface of the cell and in the IS 
Surface(  double density , double width , string cortical_Sliding , double density_in_IS , ISCorticalSl ISCorticalSl_arg );
Surface(  double density , double width , double density_in_IS , ISCorticalSl ISCorticalSl_arg , double density_in_IS_2 , ISCorticalSl ISCorticalSl_arg_2 );






//Copy constructor
Surface(const Surface& other);
//Overloading of the operator
Surface& operator=( const Surface &tmp );



//Returns the size of one bin in the map along x-Axis
double get_X_width(){ return this->x_width; };
//Returns the size of one bin in the map along y-Axis
double get_Y_width(){ return this->y_width; };
//Returns the size of one bin in the map along z-Axis
double get_Z_width(){ return this->z_width; };
//Returns density
double get_density(){ return this->density; };
//Returns the numbers of dyneins originally placed into capture-shrinkage IS 
double get_original_number_capt(){ return this->original_number_capt; };
//Returns the numbers of dyneins originally placed into cortical sliding IS1 
double get_original_number_cort_1(){ return this->original_number_cort_1; };
//Returns the numbers of dyneins originally placed into cortical sliding IS2
double get_original_number_cort_2(){ return this->original_number_cort_2; };
//Returns the key of the compartment containing the point projected on the plasma membrane
unsigned int get_dynein_compartment_id_projected( Vector3d point_position );
//Dimension is the constant necessary to get the right key of the compartment
unsigned int get_neccessary_dimension(){ return this->neccessary_dimension; };
//Returns the map containing the points
std::map< unsigned int , std::vector<Vector3d> > get_map(){ return this->surface; };
//Inserts the vector of points in the compartment with the key
void set_dynein_points( unsigned int key , std::vector<Vector3d> vectors_arg );
//Gets all dynein motors from the map
std::vector<Vector3d> get_all_dynein_point();
//Get all dyneins in the compartment with the key
std::vector<Vector3d> get_dynein_points( unsigned int key );
//Puts another point in the compartment
void add_dynein_point( unsigned int key , Vector3d vector_arg );
//Gets the number of dynein in the map
unsigned int get_dynein_motors_number();



//It controls, whether the point is contained in the map
bool dynein_surface_catch_control( Vector3d point_position , Vector3d& Catching_point  );
//Erases the point in the compartment with the map_key and array_index
void erase_dynein_point( unsigned int map_key , unsigned int array_index );
//Erase all points in the compartment with the key
void erase_vector_points_with_key( unsigned int map_key );
//Creates the dyneins for the capture-shrinkage IS
//All the points from the map are erased
void erase_map( );
std::vector<Vector3d> create_capture_shrinkage_points(  unsigned int number_of_points  , IS_Capt_Shrinkage IS_tmp );
//Adds points to the compartment specified by compartment_id
void add_dynein_points_to_compartment_with_number( std::vector<Vector3d> vectors , unsigned int compartment_id );
//Assures that the point is located on plasma membrane
Vector3d project_point_on_surface( Vector3d position );
//Adds points to the surface
void project_and_add_points_to_surface( std::vector< Vector3d > points );
//Creates cortical sliding dyneins
std::vector<Vector3d> create_cortical_sliding_points(  unsigned int number_of_points  , ISCorticalSl IS_tmp );

//This function assures the uniform distribution of the dyneins in the cortical sliding IS
void change_part_of_IS_points( ISCorticalSl IS_tmp );
//This function assures the uniform distribution of the dyneins in the cortical sliding IS
void change_part_of_two_IS_points(  ISCorticalSl IS_tmp_1 , ISCorticalSl IS_tmp_2 );
//This function assures the uniform distribution of all dyneins on the membrane
void change_part_of_points_surface_and_IS( ISCorticalSl IS_tmp );
//Generate the number from the triangular distribution between zero and upper boundary
double triangular_distribution( double upper_boundary );
//It prints the map to the terminal-just for very basic checking
void print_inner_map();

//Gets the number of dynein points in both IS - just for control
std::vector<unsigned int> get_number_of_point_in_both_IS(  ISCorticalSl IS_tmp_1 , ISCorticalSl IS_tmp_2 );


~Surface();

//////////////////////////////////////////////////////////
//Constructors for cortical sliding IS just for the purposes of testing
Surface(  double density , double width , string cortical_Sliding , double density_in_IS );
Surface( string cortical_Sliding , ISCorticalSl ISCorticalSl_arg , unsigned int minus_integer );
//This is just a test of intersection, not used for simulations producing numerical results
bool dynein_surface_catch_control_tangent( Vector3d point_position , Vector3d a_point , Vector3d tangent_b , Vector3d& Catching_point );
//This function generates dynein motors outside IS - so far, not used to get numerical results
std::vector<Vector3d> create_Cortical_Sliding_points_outside_IS( unsigned int number_of_points , ISCorticalSl ISCorticalSl_arg );




private:
    double A_axis;
    double B_axis;

    double x_width;
    double y_width;
    double z_width;
    unsigned int original_number_capt;
    unsigned int original_number_cort_1;   
    unsigned int original_number_cort_2;         
    
    const unsigned int neccessary_dimension = 1000;

    double density; //number of dynein per micro^2
    std::map< unsigned int , std::vector<Vector3d> > surface;

};

#endif // SURFACE_H
