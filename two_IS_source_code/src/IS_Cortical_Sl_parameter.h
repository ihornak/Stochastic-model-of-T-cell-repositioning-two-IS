#ifndef IS_CORTICAL_SL_PARAMETER_H_
#define IS_CORTICAL_SL_PARAMETER_H_
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
//#include <unsupported/Eigen/IterativeSolvers>
#include <eigen3/Eigen/SparseLU>
//#include <eigen3/Eigen/src/SparseLU>
using Eigen::MatrixXd;
using namespace Eigen;

namespace IS_Cortical_Sl_parameter
{
	//First cortical sliding IS
	//The distance between the center of the cell and the center of the front side of the IS
        const double front_radius = 3.8e-6;
	//The distance between the center of the cell and the center of the rear side of the IS        
        const double rear_radius = 5.5e-6;
        //Outer radius of the IS
        const double radius = 2e-6;
        //Inner radius - it would serve to model the IS as the ring
        //So far it was not used: radius_inner = 0
        const double radius_inner = 0;

	//Second cortical sliding IS	
        const double front_radius_2 = 3.8e-6;
        const double rear_radius_2 = 5.5e-6;
        const double radius_2 = 2e-6;
        const double radius_inner_2 = 0;

/////////////////////////////////////////////////////////////////////////
	//Not used in current simulations
	const unsigned int number_of_polygon = 12;
	const unsigned int number_of_mito = 5;    
	const double force_per_lenght = 5e-6;//N/microM
	const double force = 50e-12;

        const double replace_boundary =2.5e-7;
        const double catching_radius = 6.8e-6;
        const double zero_force_boundary = 3e-7;
        const double ellipsoid_control = 0.05;
        const Vector3d orientation( 0.0 , 0.0 , - 1.0 );
    
}






#endif /* IS_CAPTURE_SHRINKAGE_PARAMETERS_H_ */
