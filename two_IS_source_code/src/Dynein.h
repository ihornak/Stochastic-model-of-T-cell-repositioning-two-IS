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

#ifndef DYNEIN_H
#define DYNEIN_H
#include <iostream>
#include <stdio.h>     
#include <math.h>
#include "simulationofCell.h"

using namespace std;

namespace Dynein
{
    const double k_d = 1.0; //s^{ -1 }
    const double F_D = 3e-12;
    
    double calculate_detach_probability( double Force );
    double calculate_detach_probability_per_time_step( double Force );
        
    
    const double F_stall = 6e-12;    
    double backward_stepping_probability_between_0_Stall_force( double Force );
    double backward_stepp_prob_smaller_than_Stall_f_per_step( double Force );//
        
    
    const double v_F = 1000.0e-9; //nm s^{ -1 }
    const double d = 8e-9;//m    
    const double forward_stepping_probability = v_F / d;
    const double forward_stepping_probability_per_step = v_F / d / ( 1.0 / sim_of_Cell::time_Step );
    
    const double v_b = 6e-9;//m * s^{-1}
    const double backward_prob_force_bigger_then_Stall_force = v_b / d;
    const double backward_prob_force_bigger_then_Stall_force_per_step = backward_prob_force_bigger_then_Stall_force / ( 1.0 / sim_of_Cell::time_Step );
    
    
    const double L_0 = 110e-9; //m
    const double tmp_cutting_distance = 200e-9;
    //const double L_0 = 2e-6; //m
    const double K_a = 5.0;
    const double attach_rate_per_step = 5.0 * sim_of_Cell::time_Step;
    const double multiplicator = 1.0;
    const double reduction_abscissa = 2e-9;
    const double L_0_perpen = 10e-9;
    
    const double alpha = 1e-4;
    
    double detachment_probability( double Force );
    double detachment_probability_per_step( double Force );
    
    const unsigned int number_of_segment_steps = 7;
    
    
    
    /*
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
