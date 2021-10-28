
#ifndef DYNEIN_REAL_H
#define DYNEIN_REAL_H
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "simulationofCell.h"

using namespace std;

namespace Dynein_real
{
    const double k_d = 1.0; //s^{ -1 }
    const double F_D = 2e-12;
    const double F_stall = 4e-12;
    const double step = 8e-9;//m
    const double v_F = 1000.0e-9; //nm s^{ -1 }
    const double forward_stepping_probability = v_F / step;
    const double forward_stepping_probability_per_step = v_F / step * sim_of_Cell::time_Step;
    const double v_b = 6e-9;//m * s^{-1}
    const double backward_prob_force_bigger_then_Stall_force = v_b / step;
    const double backward_prob_force_bigger_then_Stall_force_per_step = backward_prob_force_bigger_then_Stall_force * sim_of_Cell::time_Step;

    const double K_a = 5.0;
    const double attach_rate_per_step = K_a * sim_of_Cell::time_Step;
    const double alpha = 4e-4;
    const double L_0 = 18e-9; //m

    double calculate_attachment_probability( double distance );
    double calculate_attachment_probability_per_time_step( double distance );

    double calculate_detach_probability( double Force );
    double calculate_detach_probability_per_time_step( double Force );



    double backward_stepping_probability_between_0_Stall_force( double Force );
    double backward_stepp_prob_smaller_than_Stall_f_per_step( double Force );//








    double detachment_probability( double Force );
    double detachment_probability_per_step( double Force );




    /*
    const double multiplicator = 1.0;
    const double reduction_abscissa = 2e-9;
    const double L_0_perpen = 10e-9;



    const double k_d_0 = 5.0; //s^{-1}
    const double k_a = 5.0; //s^{-1}
    const double Stall_force = 1.2e-12;
    const double Beta = 1.0 / ( 300.0 * 1.3806503e-23 );


    //extern double k_d;
    const double ATP_concentration = 1;
    double calculate_delta( );
    double calculate_k_d( double Force );

    const double v_b = 6e-9; //m/s
    const double d = 8e-9; //m
    const double backward_rate = v_b / d;

    //double calculate_Stall_force();
    */
}

#endif // DYNEIN_H
