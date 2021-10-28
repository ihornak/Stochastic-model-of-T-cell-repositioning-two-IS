#include <iostream>
#include "src/MTOC2.h"
#include "src/Microtubule.h"
#include "src/Cell.h"
#include "src/Cell_two_IS.h"
#include "ISCorticalSl.h"
#include "generalUtilities.h"
#include "Surface.h"
#include "Dynein.h"
#include "numerical_result.h"
#include "IS_Cortical_Sl_parameter.h"
#include "Generator.h"






int main(int argc, char **argv )
{
     //////////////////////////////////////////////////////////////////////////////////////////////// 
     //This part performs parallel computation of multiple simulation runs and saves numerical results
     cout.precision(17);
    //This is the initialization of the random number generator located in mersenne.h/cpp files
    //We use Mersenne Twister and the generator is initialized with time
    //Every thread has its onwn generator, we use 28 threads(increase if the number of threads > 28)
    unsigned int number_of_threads = 28;
    mersenne_twisters_64::initialize_generators( number_of_threads );
   
   

    //////////////////////////////////////////////////////////////////////////////////////////////// 
    //This part produces the files for the video depicting the process
    //The number of microtubules
    unsigned int number_Of_Microtubules = sim_of_Cell::Microtubules;
    //The cytoskeleton is symetrical - no extra microtubules are created at one side of the MTOC
    unsigned int number_Of_extra_Microtubules = 0;
    //The duration of the simulation 
    double Time_0 = 0.0;
    double Time_1 = sim_of_Cell::Time_1;
    //The number of pictures for the video of the repositioning
    unsigned int number_of_pictures = 10;
    //This is an identifier of the process - useless in the case of only one simulation run
    unsigned int index = 0;

    //Constructor of the cell
    //Creation of the cell with two IS 
    //The cytosleton is not in the configuration of the lowest energy.
    Cell_two_IS cell( number_Of_Microtubules , number_Of_extra_Microtubules , index );
    //Therefore, 10 seconds are given to find physically realistic confifguration
    //Microtubules remain free, unatached to dynein
    cell.timeDevelopment_two_Cell_calib( Time_0 , 10 ); //Time_1

    //One simulation run
    cell.timeDevelopment_Cell_two_IS( Time_0 , Time_1 , number_of_pictures );








    //This part of the code generates numerical results
/*  
    //measurement of the time
    time_t start,end;
    time ( &start );
   
    //Parallel computation of the process 
    //the details can be found in the file numerical_result.cpp
    numerical_result::test_parallel_two_IS( );
   
    //measurement of the time and printing of the duration of the simulation
    time(&end);
    double dif = difftime( end , start );
    printf ("Elasped time is %10.6lf seconds.\n", dif );
*/     
    
    
    
    
    
 
 
 
   
  
    return 0;

}


