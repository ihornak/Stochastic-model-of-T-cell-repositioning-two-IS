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

#include "Dynein.h"

namespace Dynein
{
    //double k_d = 5;
    
    double calculate_detach_probability( double Force )
    {
        return k_d * exp ( fabs( Force ) / F_D );
    }
    
    double calculate_detach_probability_per_time_step( double Force )
    {
        double probability_per_second = calculate_detach_probability( Force );
        double constant_per_time_step = ( 1.0 / sim_of_Cell::time_Step );
        return probability_per_second / constant_per_time_step;
    }
    
    double backward_stepping_probability_between_0_Stall_force( double Force )
    {
        /*
        cout<<"Force / F_stall = "<<Force / F_stall<<endl;
        cout<<"v_F / d = "<<v_F / d<<endl;
        cout<<"v_F / d * ( 1.0 - Force / F_stall ) = "<<v_F / d * ( 1.0 - Force / F_stall )<<endl;
        */
        return v_F / d * ( 1.0 - Force / F_stall );        
    }
    
    double backward_stepp_prob_smaller_than_Stall_f_per_step( double Force )
    {
        double probability_per_second = backward_stepping_probability_between_0_Stall_force( Force );
        double constant_per_time_step = ( 1.0 / sim_of_Cell::time_Step );
        return probability_per_second / constant_per_time_step;
    }
    

    double detachment_probability( double Force )
    {
        return k_d * exp( fabs( Force ) / F_D );
        //return 0.0;
    }
    
    double detachment_probability_per_step( double Force )
    {
        //cout<<"Force = "<<Force<<endl;
        double prob = detachment_probability( Force ) / ( 1.0 / sim_of_Cell::time_Step );
        prob = 1;
        //cout<<"prob = "<<prob<<endl;
        return prob;    
    }
    
    
    
    
/*    
    double calculate_delta(  )
    {
        double argument_1 = 2.9e5 / pow ( ATP_concentration , 0.3 ) - 28111.1;
        double argument_2 = 8534.6;
        if( argument_1 > argument_2 )
        {
            return argument_1 * 1e-9;
        }
        else if( argument_1 > argument_2 )
        {
            return argument_2 * 1e-9;
        }        
    }
    
    double calculate_k_d( double Force )
    {
        cout<<"Beta = "<<Beta<<endl;
        cout<<"Force = "<<Force<<endl;
        cout<<"calculate_delta( Force ) = "<<calculate_delta( )<<endl;
        double exponent = Beta * Force * calculate_delta( );
        cout<<"exponent = "<<exponent<<endl;
        cout<<"exponent = "<<exp( exponent )<<endl;
        return exponent;
    }
  
*/    
    
    
}
