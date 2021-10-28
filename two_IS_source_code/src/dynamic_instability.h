#ifndef DYNAMIC_INSTABILITY_H_
#define DYNAMIC_INSTABILITY_H_
#include "Cell_parametres.h"
#include "simulationofCell.h"


//In this header all the parameters of the dynamic instability can be found


namespace dynamic_instability
{
    //Growing speed
    const double polymerization_constant = 0.1e-6;
    //Shrinking speed
    const double shrinking_constant = 0.2e-6;
    //Rescue rate
    const double rescue_rate = 0.044;
    //Basic catastrophe rate
    const double catastrophe_rate = 0.022;    
    //Parameters describing lenght dependent catastrophe rate 
    const double coefficient_5 =  sim_of_Cell::PI * Cell_parametres::A_AXIS;
    const double coefficient_6 = sim_of_Cell::PI * Cell_parametres::A_AXIS + Cell_parametres::A_AXIS / 2.0;    
    //Computes the catastrophe rate depending on the length
    double rate_based_on_Length_August( double L_tmp );
    //These parameters are used in the function simulating the effects of the dynamic instability to obtain microtubule lengths corresponding to the correct probability distribution
    const double  Time_Of_initial_process = 600;
    const double Time_Step_Of_in_pr = 0.01;	
    //If false, normal process of repolarization
    //If true, probability density is created 
    const bool calibration_switch = false;
    //Decides whether dynamic instability is implemented in the construtor
    const bool instability_switch = true;
    //Decides whether the dynamic instability is implemented
    const bool instability_switch_2 = true;
    //Only the last segment of the microtubule can attach to capture-shrinkage dynein 
    const bool capture_shrinkage_switch = true;	    
    //Parameter chosing the dynamic instability described in the publication
    const unsigned int decision = 2;	

   //Paramters for different "dynamic instabilities"
    const double coef_leng_2 = 0.5e-6;		
    double rate_based_on_L( double L_tmp );
    const double coefficient_3 = 1.082254;
    const double coefficient_4 = 0.4e-6;
    const double coefficient_1 = 0.8;
    const double coefficient_2 = 1.8;    
    const double L_zero = coefficient_1 * sim_of_Cell::PI * Cell_parametres::A_AXIS;
    const double L_catastrophe = coefficient_2 *  sim_of_Cell::PI * Cell_parametres::A_AXIS;    

}

#endif





