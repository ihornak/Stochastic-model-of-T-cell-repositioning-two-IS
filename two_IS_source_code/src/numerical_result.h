#ifndef NUMERICAL_RESULT_H
#define NUMERICAL_RESULT_H


#include <iostream>
#include <random>

#include "MTOC2.h"
#include "Microtubule.h"
#include "Cell.h"
#include "Cell_two_IS.h"
#include "ISCorticalSl.h"
#include "generalUtilities.h"
#include "Surface.h"
#include "Dynein.h"
#include "IS_Capture_shrinkage_parameters.h"
#include "MTOCparam.h"
#include "mersenne.h"
#include <random>
#include <iostream>





namespace numerical_result
{
    //This is the simulation of one process
    void repolarization_sample_two_IS( unsigned int index );
    //This function only executes all the simulation runs given by repolarization_sample_two_IS
    void test_parallel_two_IS( );
    //The same for the cell with one IS 
    void create_document_file();
    //Creates the file with the parameters of the simulation
    void  create_document_file_two_IS();
}






#endif // NUMERICAL_RESULT_H
