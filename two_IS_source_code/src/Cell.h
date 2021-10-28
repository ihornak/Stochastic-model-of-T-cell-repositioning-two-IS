/*
 * Cell.h
 *
 *  Created on: Mar 10, 2016
 *      Author: hornak
 */

#ifndef CELL_H_
#define CELL_H_

#include "simulationofCell.h"
#include "Microtubule.h"
#include "Nucleus.h"
#include "MTOC2.h"
#include <errno.h>
using Eigen::MatrixXd;

#include <omp.h>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include "ISCaptShrinkage.h"
#include "IS_Capture_shrinkage_parameters.h"
#include "ISCorticalSl.h"
#include "MTOCparam.h"
#include "CellShape.h"
#include "Cell_parametres.h"
#include "IS_Cortical_Sl_parameter.h"
#include "IS_Dynein_Cell_surface.h"
#include "GeneralGeometry.h"
#include "Surface.h"
#include "mtoc_simple.h"
#include "ideal_mtoc.h"
#include <fstream>               //mandatory for ofstreaming
#include <sstream>
#include <random>
#include <iostream>

using namespace std;



class Cell
{
public:
	//This is constructor for unconstrained cytoskeleton in 3D space
    	//Deafault constructor	
    	Cell(  );
    	Cell( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant );


     //Getters Setters	
    	//Returns the microtubule identified by the index
    	Microtubule getMicrotubule( unsigned int index );
    	//Returns the number of microtubules in the cell    
    	unsigned int get_microtubule_number();
    	//Returns the capture-shrinkage
    	IS_Capt_Shrinkage get_IS_capture_shrinkage();  
    	//Estimates the friction of the entire cytoskeleton
    	double get_cytoskeleton_friction();
    	//Gets the number of attached capture-shrinkage dyneins
    	unsigned int get_number_of_dynein_motors_in_micro_capture();
    	//Returns cortical sliding IS
    	ISCorticalSl get_IS_cortical_sliding_first(){ return this->IS_Cortic_Sl; };
    	//Sets the map of the capture-shrinkage IS
    	void set_density_IS_capture_shrinkage( Surface& density ){ this->Capture_Shrinkage_dynein = density; };    	
    	//Sets the map holding cortical sliding dyneins
    	void set_density_surface_dynein( Surface& density );       	
         //Get number of cortical sliding dyneins attached to microtubules
        unsigned int get_number_of_dynein_motors_cortical_sliding_micro();   	

  	//MICRO WALL NUCLEUS INTERACTION
        //Simplied version of the force, the cell is taken as a sphere
        Vector3d force_position_cell_wall_sphere( Vector3d position  );    	
    	//The force of the wal acting on the microtubule bead
        Vector3d force_position_cell_wall_two_elipsoid( Vector3d position  );
        //Computes the force acting on one point of the MTOC - force_position_cell_wall_two_elipsoid
        Vector3d MTOC_point_wall_interaction( Vector3d position );           
	//The force from the wall acting on every bead of the microtubule	
        void force_whole_microtubule_cell_wall( MatrixXd &force , unsigned int index );
        //Computes the force of the wall on  the MTOC
        void force_MTOC_cell_wall( MatrixXd &force );         
        //Computes the force of the nucleus acting on the entire microtubule given by index
        void force_whole_microtubule_nucleus( MatrixXd &force , unsigned int index ); 
        //Force of the MTOC acting on all MTOC points
        void MTOC_nucleus_interaction( MatrixXd &force_on_MTOC );    
        
     //Micro-MTOC interaction    
        //Interaction between microtubule and the MTOC
        void MTOC_microtubule_two_points_force_with_bending( Vector3d& first_point_force , Vector3d& second_point_force , Vector3d& bending_opposite_MTOC_point , Vector3d& bending_MTOC_point , unsigned int microtubule_number );
	//Simplified calculation of the force between the microtubule and MTOC
        Vector3d MTOC_microtubule_two_points(  unsigned int microtubule_number );

        
    //MICRO DYNEIN CELL IS INTERACTION
 	//Returns the all dyneins in the compartment determined by ID_map_index
        std::vector<Vector3d> get_dynein_in_compartment( unsigned int ID_map_index );
       //Cortical sliding dynein attaches on microtubule
        void catch_pair_abscissa_real_dynein();
	//Adds new dynein to microtubule
        void microtubule_catch_pair_abscissa_real_dynein( unsigned int microtubule );
        //This function resizes the microtubule and performs the stepping of dynein
        void stepping_and_detachment_of_all_microtubule_projection_real_dynein();
        //This is the function for control, whether the microtubule 0 or 9 is caught
        void check_caught_micro_IS_with_real_dynein( unsigned int microtubule );
	//Checks whether a new microtubule is attached to capture-shrinkage dynein
        void check_caught_micro_IS_with_real_dynein_all_micro( );
	//Adds new capture-shrinkage dyneins
	void microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro( );
        //Adds new dynein to the microtubule already attached to the first IS
        void microtubule_catch_pair_abscissa_real_dynein_in_IS( unsigned int microtubule );
        //Checks the length of microtubule because of the numerical inprecisions during attachment
        void control_length_of_micro_IS();
        //checks the length of attached microtubules
	void control_length_of_micro_IS_2();
////////////////////////////////////////////////////////////////////////////////////






   //general utilities
       //It takes the point and projects it on plasma membrane
       Vector3d project_point_on_surface( Vector3d position );
       //It returns the id of the compartment containing the point
       unsigned int get_dynein_compartment_id( Vector3d position );
       //Projects the points on plasma membrane and adds them to the map
       void project_and_add_points_to_surface( std::vector< Vector3d > points );
       //Adds points to the map
       void add_points_to_surface( std::vector< Vector3d > points );
	//Resizing of the cytoskeleon due to numerical imprecisions
        void resize_micro_with_different_first_segment_and_MTOC(  );






    	//PRINTING
    	//Saves the configuration of the cell for the purposes of visualisation
        void print_Cell( unsigned int index );
        //Saves the basic shapes of the cell
        void print_Cell_parametres( );
        //Saves the time and the center of the MTOC            		
        void BIG_PRINT( unsigned int time_Step );
        //Saves the center of the MTOC, number of attached capture-shrinkage and cortical sliding dyneins into files with the name containing the id of the process
        void BIG_PRINT_numerical( unsigned int time_Step , unsigned int key );
        //Saves the positions of all microtubule beads and all dyneins
        void BIG_PRINT_numerical_micro( unsigned int time_Step , unsigned int key );
        //Saves the time into the file 
        void print_time( unsigned int time_Step );
        //Saves the center of the MTOC into the file
        void print_center_MTOC();
	//Saves the center of the MTOC into the file containing the id key of the process
        void print_center_MTOC_numerical( double time_Step , unsigned int key );
        //Saves the number of attached capture-shrinkage dyneins into a file containing the name of the process
	void print_dynein_capture_shrinkage_numerical( double time_Step , unsigned int key );	
	//Saves the number of cortical sliding dyneins acting on MTs
        void print_dynein_cortical_sliding_numerical( double time_Step , unsigned int key );           
        //This function saves the entire configuration of all microtubules in cytoskeleton
	void print_microtubules_numerical( unsigned int time_Step , unsigned int key );
	//Saves the positions of all dyneins - attached ot unattached
 	void print_dynein_numerical( double time_Step , unsigned int key );
	//This saves the front center of the IS
	void print_center_numerical_parameters( unsigned int key ); 	

	//The Mid-Step algorithm
        void MidStep_3();







    //void one_simple_time_step();

    void timeDevelopment_2( double Time_0 , double Time_1 );
    void timeDevelopment_2( double Time_0 , double Time_1 , unsigned int number_of_Pictures );


    void timeDevelopment_numerical_results_2( double Time_0 , double Time_1, unsigned int key );
    
    //Computes the process with the Euler algorithm - unused, just for cheking
    void Euler_algorithm( );
    void Euler_timeDevelopment( double Time_0 , double Time_1 );
    void Euler_timeDevelopment( double Time_0 , double Time_1 , unsigned int number_of_Pictures );    
    


	virtual ~Cell();
protected:
	double a_axis;
	double b_axis;
	Nucleus nucleus;
        MTOC2 MTOC;

        unsigned int number_of_microtubules;
        unsigned int number_of_microtubules_extra;
        Microtubule *array_Of_Microtubules = NULL;
        IS_Capt_Shrinkage IS_Capture_Shrinkage;
        ISCorticalSl IS_Cortic_Sl;
        Surface density_surface_dynein;
        Surface Capture_Shrinkage_dynein;

        ideal_MTOC abstract_MTOC;
        MTOC_simple mtoc;
        std::random_device rd;


};

#endif /* UNCONSTRAINEDCELL_H_ */
