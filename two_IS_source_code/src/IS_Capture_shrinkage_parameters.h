/*
 * IS_Capture_shrinkage_parameters.h
 *
 *  Created on: Apr 7, 2017
 *      Author: hornak
 */

#ifndef IS_CAPTURE_SHRINKAGE_PARAMETERS_H_
#define IS_CAPTURE_SHRINKAGE_PARAMETERS_H_
#include "simulationofCell.h"


//IS_Capture_shrinkage_param::depolimerization_treshold
namespace IS_Capture_shrinkage_param
{
   //Parameters describing the positions and sizes of both capture-shrinkage IS
////////////////////////////////////////////////////////////////////////////////////////////////
    //First capture-shrinkage
    //Radius of the IS
    const double radius = 0.4e-6;
    //Azimuthal angle of the IS
    const double azimutal_angle = 0.0;
    //Polar angle of the IS
    const double polar_angle = 3.0 / 8.0 *  sim_of_Cell::PI;   //3.5 / 8.0
    //Distance between the center of the cell and the center of the front side of the IS
    const double z_coordinate_front = 3.7e-6;
    //Distance between the center of the cell and the center of the rear side of the IS
    const double z_coordinate_back = 6.3e-6;
////////////////////////////////////////////////////////////////////////////////////////////////
    //Second  capture-shrinkage
    const double radius_2 = 0.4e-6;
    const double azimutal_angle_2 = 1.0 *  sim_of_Cell::PI;
    const double polar_angle_2 = 3.0 / 8.0 *  sim_of_Cell::PI;   //3.5 / 8.02
    const double z_coordinate_front_2 = 3.7e-6;
    const double z_coordinate_back_2 = 6.3e-6;
///////////////////////////////////////////////////////////////////////////////////////////////





    //Unused	
    const unsigned int number_of_polygon_lower = 0;
    const unsigned int number_of_polygon_upper = 1;
    const unsigned int number_of_mito = 1;
    const double cut_distance = 0.5e-6;
    const double dynein_Force_capture_shrinkage = 4e-12;
    const double cutting_distance_depol_micro = 0.3e-6; 
    const double depolimerization_distance_threshold = 2e-6;
    const double procentage_constant = 1.3;
    const double minimal_radius =  4.9999e-6;








    const double depolimerization_treshold = 1e-7;


}






#endif /* IS_CAPTURE_SHRINKAGE_PARAMETERS_H_ */
