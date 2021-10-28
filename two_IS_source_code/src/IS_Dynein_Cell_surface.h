#ifndef IS_DYNEIN_CELL_SURFACE_H_
#define IS_DYNEIN_CELL_SURFACE_H_

//This would be used if the dynein was distributed on the cell mebrane outside IS
//It is not used in current simulations
namespace IS_Dynein_Cell_surface
{
	const unsigned int number_of_polygon_lower = 0;
    	const unsigned int number_of_polygon_higher = 12;
	const unsigned int number_of_mito = 5;    
	const double force_per_lenght = 0.18e-6;//N/microM 0.16   0.12e-6
	const double force = 5e-12;//N/microM 0.45  10e-12
	const unsigned int lower_side = 2;
    	const double close_boundary = 0.3e-6;
        const double dynein_on_surface_catch_radius = 12.5e-9;
        const double treshold_cut_tangent = 10.0e-7;
         const double new_radius = 0.075e-7;
    
}
//IS_Dynein_Cell_surface::force

#endif
