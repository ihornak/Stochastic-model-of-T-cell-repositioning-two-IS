/*
 * microtubule.h
 *
 *  Created on: Jul 24, 2015
 *      Author: hornak
 */

#ifndef MICROTUBULE_H_
#define MICROTUBULE_H_
#include "MTOCparam.h"
#include "KISS.h"
#include <fstream>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <map>
#include "mersenne.h"
#include <random>
#include <iostream>


#include "simulationofCell.h"
#include "ISCaptShrinkage.h"
#include "Cell_parametres.h"
#include "IS_Cortical_Sl_parameter.h"
#include "IS_Capture_shrinkage_parameters.h"
#include "GeneralGeometry.h"
#include "IS_Dynein_Cell_surface.h"
#include "Dynein.h"
#include "Dynein_real.h"
#include "dynamic_instability.h"



#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>	
#include <Eigen/SparseQR>
#include <eigen3/Eigen/IterativeLinearSolvers>
using Eigen::MatrixXd;
using namespace Eigen;
typedef Eigen::Triplet<double> T;

using namespace std;



class Microtubule {
public:

    	//This constructor is used in the simulations
        Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points  );       
 	//Microtubules in unconstrained 3D space
	Microtubule();
   	//Copy constructor
	Microtubule( const Microtubule & tmp );



    //SETTERS GETTERS
        //Get number of mirotubule beads
	unsigned int getNumberOfPoints()  const;
	//Sets the number of beads of the microtubule 
	void setNumberOfPoints( unsigned int number_of_point_tmp ){ this->numberOfPoints = number_of_point_tmp; };
	//Sets the coordinates of the microtubule
        void setCoordinates( MatrixXd &coordinatesTmp );
        //Returns the coordinates of the microtubule
        MatrixXd getCoordinates() const;    
        //Sets the position of one bead    
        void setPoint( unsigned int number_of_point ,  Vector3d position );
	//Sets the length of the segment
	void setRestDist( double tmp_rest ){  this->restDistancePoints = tmp_rest; };
	//Returns the length of the segment
        double getRestDist() const;
        //Returns the bending rigidity of the discretized polymer
        double getKappa() const;
        //Sets the friction of the bead
	void set_effective_friction( double tmp ){ this->effective_friction = tmp; };
	//Returns the effective friction of the microtubule bead
        double get_effective_friction();
        //Returns the effective friction of the whole microtubule
        double get_effective_friction_whole_microtubule();
        //Returns the position, where the microtubule depolimerizes
        Vector3d get_IS_position_catching( ) const;
         //Returns the identifier of the microtubule
        unsigned int getID()  const;
        //Returns the side of the MTOC where the microtubule is anchored
        unsigned int getSide()  const;
        //Returns the number of MTOC polygon, to which the microtubule is attached 
        unsigned int get_polygon_number()  const;
        //Returns the MTOC point from which the microtubule is sprouting 
        unsigned int get_MTOC_point();
        //Return rear points of the MTOC
        unsigned int get_MTOC_opposite_point();
        //Determines whether the microtubule is free, or attach to capture-shrinkage or cortical sliding dyneins
        unsigned int get_dynein_index() const;
	//It does what it says
        void set_dynein_index( unsigned int new_value );
        //Sets the point, where the microtubule is depolimerized
        void set_IS_position_catching( Vector3d IS_vector );
        //Get the point determined by the index
        Vector3d getPoint( unsigned int index )  const;
        //Returns the last point of the microtubule
        Vector3d get_last_Point( )  const;
        //Returns the segment connecting two beads
        Vector3d getTangent2( unsigned int index ) const;
        //Gets the tangent connecting the last two points of the microtubule
	Vector3d get_last_Tangent( )  const;
	//Prolongs the last segment of the microtubule by additional_lenght
	void prolong_Last_tangent( double additional_lenght );
	//Shortens the last segment of the microtubule
	void shorten_Last_tangent( double shortening );
	//Gets the length of the microtubule
        double get_lenght_of_microtubule() const; 
        //Gets the length of the microtubule outside of MTOC
	double get_lenght_of_microtubule_outside_MTOC() const;
	//The lenght of the microtubule till the bead with the index 
        double get_distance_to_lower_bead_with_index( unsigned int index ) const;
        //Returns the lengths of all tangents
        MatrixXd get_lenght_of_tangents(){ return this->lenght_of_tangents; };
        //Sets the length of the microtubule after attachment to capture-shrinkage dynein
        void set_lenght_after_catching( double lenght_after_catching_arg );
        //Returns the length that the microtubule had after attachment to capture-shrinkage dyneins
        double get_lenght_after_catching( );
         //overloading of the operator
	Microtubule& operator=( const Microtubule &tmp );
	//Calculate effective friction  using mentioned publication  David Leith - Aerosol Science and Technology
	double calculateEffectiveFriction( );
	//Caculates effective friction Howard, J., 2001. Mechanics of Motor Proteins and the Cytoskeleton
	double calculateEffectiveFriction_Howard( );
	//Retruns the number of points attached to the microbule
        std::vector<Vector3d> get_dynein_points_Cortical_Sliding();	
        //Alters the grow index of the microtubule
	void set_growth_index( unsigned int index_tmp ){ this->growth_index = index_tmp; };
	//Gets the grow index of the microtubule
	unsigned int get_growth_index( ){ return this->growth_index; };
	//The microtubule lenght outside the MTOC
	//Sets the length outside of MTOC
	void set_lenght_of_micro_0_outside_MTOC( double tmp_length ){ this->lenght_of_micro_0_outside_MTOC = tmp_length; };
	//Returns the length outside of MTOC
	double get_lenght_of_micro_0_outside_MTOC(  ){ return this->lenght_of_micro_0_outside_MTOC; };
	//Saves all the lengths of all segments connecting all the beads
	void set_lenght_of_tangents_micro_0(){ this->lenght_of_tangents_micro_0 = lenght_of_tangents; };
        //Simply removes the last bead of the microtubule
    	void set_Cell_ID( unsigned int id_tmp ){ this->cell_id = id_tmp; };
       	//Returns the number of points attached to the microbule
        unsigned int get_number_of_dynein_points_IS();    	
  
     //Printing of microtubule
	//This prints the beads into terminal - basic control
        void print_Microtubule();
         //This prints the tangents of microtubules into terminal - basic control
        void print_Microtubule_Tangent();    	
        
     //Projection
	//Creates the matrixes for the projection of forces
	void getMatrixes( MatrixXd &inv_Matrix , MatrixXd &projection_Matrix );
	//Creates the matrix expresing derivatives of constraints
	void getN_i_mu( MatrixXd &n_i_mu );
	//Gets the projection matrix and projects forces
	void project_force( MatrixXd& original_force , MatrixXd& projected_forces );         
    	
    	//Simply removes the last bead of the microtubule
    	void Substract_one_Bead_simple();
    	// It returns random forces acting on the whole microtubule
	void getRandomForces(  double timeStep , MatrixXd &randomForces );
	//Computes bending forces
	void getBendingForces( MatrixXd &bendingForce );
	//Bending force calculation using matrixes
        void getBendingForces_2( MatrixXd &bendingForce );
        //It determines whether the point is in the cell 
        bool confirm_inner_ellipsoid( Vector3d position , double a_axis , double b_axis );
        //Analytical calcul
 	//Function deletes the last bead of microtubule
	void Substract_one_Bead();   
	//Controls whether NaN exists in the variables of the microtubule             
        void control_NAN(  string arg );
        //Controls whether NaN exist in the matrix
	bool control_nan_Utilities( MatrixXd tmp );
	//Places the last bead of the microtubule to the position, where microtubule depolymerizes
        void add_one_final_Bead_catching_position();
        //The force acting on a point on the segment is transmitted to two beads at the end of segments
        void distribute_force_on_beads( MatrixXd& force_on_microtubule , Vector3d force , unsigned int bead_segment , double t );

      //Dynein forces
	//Returns the simplified force of dynein - costant force
	Vector3d dynein_force_capture_shrinkage( );
        //This adds one pair - anchor point and abscissa to microtubule    
        void add_pair( std::pair < Vector3d ,double  > point_abscissa );
        //Returns all the pairs(anchor point, abscissa) of dynein acting on microtubule  
        std::vector< std::pair < Vector3d , double > > get_Dynein_abscissa(){ return this->Dynein_motors_2; }; 
        //Returns all dyneins acting on microtubules and erase them       
        std::vector< Vector3d  > get_Dynein_points_and_erase(  );
        //The same - intentional duplication
        std::vector<Vector3d> get_dynein_points_in_IS_and_erase();
        //Returns all the dyneins acting on microtubule
        std::vector< Vector3d  > get_Dynein_points_without_erasing(  );  
        //Matrix contain all the forces from the dynein
        void force_dynein_abscissa_real_dynein( MatrixXd& force_dynein );      
         //Returns force of dynein acting on microtubule        
        Vector3d force_real_dynein_one_pair(  std::pair < Vector3d , double > pair_tmp   );
        //This function performs the steppping and detachment of dynein on the microtubule
        std::vector< Vector3d > stepping_detach_real_dynein_abscissa_projection( ); 
        //Returns the number of segment(or lower bead), in which abscissa is located
        unsigned int get_index_according_to_abscissa( double abscissa );
        //Returns attachment point according to the abscissa
        Vector3d get_attachment_point_according_to_abscissa( double abscissa ); 
        //The lenght of all tangent are saved
        void set_lenght_of_tangents( string argument );
        //Controls the length of the microtubule and compares it with the distance between the MTOC and the depolymerization point 
	std::vector<Vector3d> control_length_of_micro_IS_2();
	//Returns all dyneins acting on a microtubule attached to cortical sliding
	std::vector<Vector3d> get_catching_points_IS_sliding( );     
	//Detaches from dyneins and the index is set to zero
	std::vector<Vector3d> detach_from_IS_2();
	//Check whether dyneins are not detached due to the depolymerization
        std::vector< Vector3d > detachmet_due_to_shrinkage_of_MT_9( );
  
	    
     //Resizing
        //This is the simplest resizing of the microtubule - all segments are set to this->restDistancePoints 
        void resizeMicrotubule();
        //Resizing, the first segment has the orientation given by the vector
        void resizeMicrotubule( Vector3d orientation );
        //Simple resizing of the microtubule with the different length of the first segment
        void resizeMicrotubule_first_segment_different( );
        //This resized the microtubule with the lengths given by the vector
        void resizeMicrotubule_with_different_tangent_lenghts( );
        //Resizes the microtubule so that all segments have the same size and the number of beads stays the same
	void resizeMicrotubule_same_number_of_beads();
        //Resizes the microtubule so that all segments have the same size, number of beads is increased
	void resizeMicrotubule_plus_one();
        //Resizes the microtubule so that all segments have the same size, number of beads is reduced
	void resizeMicrotubule_minus_one();        
        //Controls whether the microtubule was overextended due to the repeated attachment 
	bool control_of_captured_MT_lenght_against_prolongation(  ); 	
  	//Resizes growing, unattached microtubule
	void control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0();
        //Shortens the microtubule if it became overextended due to repeated attachment
	void shorten_MT_due_to_numerical_prolongation_during_attachment(  );	
        //Control of the microtubule length and the stepping of dynein  
        std::vector<Vector3d>control_and_resizing_and_motor_adjustment_depolimerizing_microtubule();
        //Checks the last segment of the microtubule and resizes it
	std::vector<Vector3d>  control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2();
	//Controls the last segment, resizes microtubule and checks the detachment of dyneins
	std::vector<Vector3d>  control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple();
	//Resizing of the microtubule for the microtubule attached to capture-shrinkage dynein
	void  control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple_after_catching();
	//Controls if the dyneins detach due to depolymerization
        std::vector<Vector3d> control_motor_detachment_with_depolimerization();       
        //Controls the the depolymerizing microtubule and resize it when necessary 
	void control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0();    
	//Controls a shrinking microtubule and resizes it
	void  control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9();
	
    //Dynamic instability
        //This is the growth of the microtubule due to the dynamic instability
        void microtubule_growth_MT_0_9();
        //Shrinking of the microtubule with a free end
        void microtubule_shrinkage_MT_0_9(  );
        //Basic polymerization of the microtubule
        void microtubule_growth_2(); 	
      
        
        //Mid step algorithm for the propagation of the microtubule
	void oneStepMidStepAlgorithm();
	//This is used when calculating movement in the cell
	//This is one step of midStep algorithm that includes external forces form MTOC and wall of the cell
	void oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces );
	//This is the first half of the mid step algorithm
	void oneStepMidStepAlgorithm_1_half_producing_random( MatrixXd& randomForces , MatrixXd& extrenal_Forces );
	//This is the second half of the mid step algorithm
	void oneStepMidStepAlgorithm_2_half( MatrixXd& original_Coordinates , MatrixXd& randomForces , MatrixXd& extrenal_Forces );


	//This function propagates the movement of the microtubule with Euler algorithm
    	void Euler_algorithm( MatrixXd& randomForces , MatrixXd& extrenal_Forces );    
    
  
	
















    //Currently unused
    
        //Unused - for the future simulations
        void set_bending_matrix( );
        void oneStepMidStepAlgorithm_2_half_projected( MatrixXd& original_Coordinates , MatrixXd& randomForces , MatrixXd& other_Forces );
	//Currently unused constructors
	Microtubule( unsigned int ID );
	Microtubule( Vector3d MTOC , unsigned int ID );
	Microtubule( Vector3d MTOC , Vector3d first_Point , unsigned int ID );
	//microtubule ID is number of microtubule in polygon
	Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID , unsigned int poly);   //PREDELAT POZDEJI!!!!!
	Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , double a_axis , double b_axis );
	Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , double a_axis , double b_axis , unsigned int number_of_points  );
        Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int number_of_points  );       
    
	Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points , bool rovna_mikrotubula );
    Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points , double first_bead_distance );
   	Microtubule( const char* filename , double rad_of_Cell , Vector3d IS , unsigned int ID );      
        
        
	virtual ~Microtubule();




private:
	// unique id of the microtubule. Its only purpose is to easje  orientation
	unsigned int numberOfPoints;    
    	unsigned int side;
    	unsigned int polygon;
	unsigned int microtubule_id;
    	unsigned int MTOC_point;
    	unsigned int second_MTOC_point;     
	double restDistancePoints;
    	double restDistancePoints_first;
	double kappa;
	double effective_friction;
	MatrixXd coordinates;
	Vector3d IS_position_catching;
	unsigned int dydein_index;
 

        std::vector< std::pair < Vector3d , double > > Dynein_motors_2;
        MatrixXd lenght_of_tangents;
        double lenght_after_catching;
	SparseMatrix<double> b_MATRIX;
	SparseMatrix<double> Kronecker_Delta;


	unsigned int growth_index;
	//tahle promena slouzi jako fixni delka mikrotubule, ktera se prodluzuje jen pri micro_growth
	//slouzi k tomu, aby se chycena MT neprodluzovala diky numerickym nepresnostem
	double lenght_of_micro_0_outside_MTOC;
	MatrixXd lenght_of_tangents_micro_0;
	unsigned int cell_id;



    
};

#endif /* MICROTUBULE_H_ */
