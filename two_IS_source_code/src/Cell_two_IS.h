/*
 * Cell.h
 *
 *  Created on: Mar 10, 2016
 *      Author: hornak
 */

#ifndef CELL_TWO_IS_H_
#define CELL_TWO_IS_H_

#include "Cell.h"
#include "numerical_result.h"
#include <iostream>

using namespace std;



class Cell_two_IS: public Cell
{
public:
    //This is constructor for an unconstrained cytoskeleton in 3D space.
    Cell_two_IS(  );
    //Constructor of the cell with the two IS
    //Number of microtubules is given by the file simulationofCell.h - 100
    //The variable number_Of_extra_Microtubules would be used to study the asymetrical cytoskeleton.
    //The results of the study were not published and number_Of_extra_Microtubules = 100    
    Cell_two_IS( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant );

  //Setters and getters
    //Basic setter, sets the Surface for the second IS capture-shrinkage    
    void set_density_IS_capture_shrinkage_second( Surface& density );
    //Gets the second capture-shrinkage IS
    IS_Capt_Shrinkage get_IS_capture_shrinkage_second(  );
    //Sets and gets the time 
    void set_time_clock( double tmp ){ this->time_clock = tmp; };
    double get_time_clock(  ){ return this->time_clock; };
    //Returns the number of capture-shrinkage dyneins in the IS 2 attached to the microtubules
    unsigned int get_number_of_dynein_motors_capture_shrinkage_2();
    //Returns the number of capture-shrinkage dyneins in the IS 1 attached to the microtubules        
    unsigned int get_number_of_dynein_motors_capture_shrinkage_1();
    //Returns the number of microtubules attached to capture-shrinkage dyneins in IS1        
    unsigned int get_number_of_microtubules_micro_capture_20();
    //Returns the number of microtubules attached to capture-shrinkage dyneins in IS2    
    unsigned int get_number_of_microtubules_micro_capture_40();
    //This function gives the number of microtubules, whose number of beads is smaller than a treshold  
    unsigned int get_number_of_short_micrtobules();
    //Returns the second cortical sliding IS
    ISCorticalSl get_IS_cortical_sliding_second(){ return this->IS_Cortic_Sl_2; };    
    //This function makes files for the generation of the video
    void print_Cell_two_IS( unsigned int index );
    //This prints the entire configuration of the cell into files that are later used to generate the videos using Povray    
    void print_Cell_two_IS_numerical( unsigned int index );
    //This function resizes the cytoskeleton after every step  in timeDevelopment_Cell_two_IS and timeDevelopment_Cell_two_IS_numerical_results
    void stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS( );
    //This function resizes the cytoskeleton after every step in the timeDevelopment_two_Cell_calib function 
    void stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS_calibration_axis_z( );

    //In  this function the new dynein attach to the microtubule previously attached in the capture-shrinkage IS
    void microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro_two_IS();
    //Adds new dynein to the microtubule already attached to the second IS    
    void microtubule_catch_pair_abscissa_real_dynein_in_IS_2( unsigned int microtubule );
    //Check whether a new microtubule attaches to capture-shrinkage dynein    
    void check_if_micro_is_caught_in_IS_two_IS( );
    //Check, whether unattached MT attaches to capture-shrinkage dynein 
    void check_caught_micro_IS_with_real_dynein_2( unsigned int microtubule );



    //Saving the results 
    //Saves results for statistical analysis
    void BIG_PRINT_numerical_two_IS( unsigned int time_Step , unsigned int key );
    //This functions saves the number of cortical sliding dyneins from both IS acting on microtubules     
    void print_dynein_cortical_sliding_numerical_two_IS( double time_Step , unsigned int key ); 
    //Saves the number of capture-shrinkage dyneins from both IS acting on microtubules    
    void print_dynein_capture_shrinkage_numerical_two_IS( double time_Step , unsigned int key );
    //Saves the number of microtubules attached in both capture-shrinkage IS     
    void print_microtubules_capture_shrinkage_numerical_two_IS( double time_Step , unsigned int key );
    //Saves the number of microtubules shorter than a treshold
    void print_number_of_shorter_microtubules_two_IS( double time_Step , unsigned int key );
    //Saves the data about unattached and cortical sliding microtubules    
    void print_MT_stats( double time_Step , unsigned int key  );	
    //Saves the lenght of all microtubules    
    void print_MT_stats_lenght( double time_Step , unsigned int key  );
    //This functions control, whether microtubules with no dynein were properly detached from the IS
    void control_dynein_detachment();
    //This function simulates the effects of the dynamic instability
    static MatrixXd dynamic_instability_MT_grow_August(  );
    //Implementation of the dynamic instability
    void dynamic_instability_August(  );
    //Controls the number of dynein in both IS
    void test_dynein( string param );     
    
    
    //Mid-step algorithm for the cell with two IS	
    void MidStep_two_IS();
    //In this function the cytoskeleton finds a configuration of the smallest energy
    //Microtubules grows and shrinks to assure that their lenghts are distributed according to specified probability density     
    void timeDevelopment_two_Cell_calib( double Time_0 , double Time_1 );
    //This function is used for one simulation run. It is used to collect basic numerical results
    //The function is used for the visualisation of the simulation
    //The argument number_of_Pictures determines, how many times the configuration of cell will be saved into text files
    //The files are later used to create Povray files, images and videos    
    void timeDevelopment_Cell_two_IS( double Time_0 , double Time_1 , unsigned int number_of_Pictures );
    //The same but producing the numerical results
    void timeDevelopment_Cell_two_IS_numerical_results( double Time_0 , double Time_1 , unsigned int key );

    //If not a number detected, the simulation is immediately stopped 	
    void control_NAN(  string arg  );

    
    //Currently unused 
    //Implementation of dynamic instability with different sets of parameters
    void all_micro_growth();
    void set_dynamic_instability();
    void dynamic_instability(  );
    MatrixXd dynamic_instability_MT_grow(  );
    void dynamic_instability_linear(  );
    MatrixXd dynamic_instability_MT_grow_linear(  );
    void control_lenght_prolongation_capture( );

    virtual ~Cell_two_IS();

private:
	unsigned int number_of_IS;
        IS_Capt_Shrinkage IS_Capture_Shrinkage_2;
	ISCorticalSl IS_Cortic_Sl_2;
        Surface Capture_Shrinkage_dynein_2; 
        double time_clock;  
	unsigned int cell_id;


};

#endif /* UNCONSTRAINEDCELL_H_ */
