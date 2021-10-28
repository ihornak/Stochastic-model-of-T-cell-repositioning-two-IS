#include "dynamic_instability.h"



namespace dynamic_instability
{

    double rate_based_on_L( double L_tmp )
    {
	double b_tmp = ( L_zero - L_catastrophe ) / log( catastrophe_rate );
	double c_rate = exp( ( L_tmp - L_catastrophe ) / b_tmp );
	double rate_return = c_rate;
	return rate_return;
    }	
    
    
    //Computes the catastrophe rate depending on the length
    //Used in simulation
    double rate_based_on_Length_August( double L_tmp )
    {
	double L_zero_tmp = coefficient_5;
	double L_catastrophe = coefficient_6;
	double b_tmp = ( L_zero_tmp - L_catastrophe ) / log( catastrophe_rate );
	double c_rate = exp( ( L_tmp - L_catastrophe ) / b_tmp );
	double rate_return = c_rate;
	return rate_return;
    }	


}

