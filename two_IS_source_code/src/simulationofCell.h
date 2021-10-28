/*
 * simulationofCell.h
 *
 *  Created on: Aug 19, 2015
 *      Author: hornak
 */
#ifndef SIMULATIONOFCELL_H_
#define SIMULATIONOFCELL_H_


namespace sim_of_Cell
{

    //////////////////////////////////////////////////////	+
        //Random force are switched of
        const bool random_force_switch = false;
        const bool random_force_switch_MTOC = false;
        const double time_Step = 1e-3;  //0.7e-3
        const double time_Step_half = time_Step / 2.0;// / 2.0 /// 2.0 //

	//Number of microtubules
        const unsigned int Microtubules = 100;
        //Time of the simulations
        const double Time_1 = 600;
        //Dynein density on the surface
        const double density_of_dynein_motor_surface = 0;

        //dynein density first IS
        const double density_of_dynein_motor_Cortical_Sliding = 0;
        const double density_of_dynein_in_IS_Capture_Shrinkage = 400;

        //dynein density second IS
        const double density_of_dynein_motor_Cortical_Sliding_2 =  0;
        const double density_of_dynein_in_IS_Capture_Shrinkage_2 = 400;
        const unsigned int number_of_samples = 500;







	//Physical constants
	const bool writing_switch = false;
	const double surface_width = 0.5e-6;
	const unsigned int MicrotubulePoints = 20;
	const double resting_distance = 0.9e-6; // 1.62e-6;  1.55e-6 11 / 2.0 0.8e-6
	const unsigned int REDUCTION = 5;
	const double chosen_prolongation_tr = 1.0e-6;
	const double k_bending =  2.3e-23 / resting_distance ;// 2.2e-23 / resting_distance 
        const double viscosity_water = 0.001002;
        const double viscosity = 0.09;
        const double deeper_viscosity = viscosity * 10.0;       
        const double k_bending_analytical =  2.3e-23 ;
        const double Temperature = 300.0;								//Kelvin
        const double boltzmann_constant = 1.3806503e-23;			//Joule / Kelvin 1.3806503e-23
        /////////////////////////////////////////////////////
        const double dynein_Force_cortical_sliding = 5e-12;
        const double PI = 3.141592653589793;
        const double angle = 0.0;
        const double multiply_friction_constant = 1;
        const bool paralelization_switch = 1;
        const unsigned int general_switch = 1;
        const bool cytoskeleton_asymetry_switch = true;
        const double surface_replacement_procentage = 0.1;

}






#endif /* SIMULATIONOFCELL_H_ */
