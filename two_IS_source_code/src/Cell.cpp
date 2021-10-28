/*
 * Cell.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: hornak
 */

#include "Cell.h"
#include <tuple>

//Deafault constructor
Cell::Cell(  )
{
	this->a_axis = Cell_parametres::A_AXIS;
	this->b_axis = Cell_parametres::B_AXIS;
        this->number_of_microtubules_extra = 0;
	this->nucleus = Nucleus();

	//initializing number of microtubules
	this->number_of_microtubules = 16;
	//create array of microtubules
	this->array_Of_Microtubules = new Microtubule[ this->number_of_microtubules ];


	this->IS_Capture_Shrinkage = IS_Capt_Shrinkage();
	//this->IS_Cortical_Sliding = ISCorticalSliding();
        this->IS_Cortic_Sl = ISCorticalSl();

        this->MTOC = MTOC2( this->number_of_microtubules );
        //this->abstract_MTOC = ideal_MTOC();
}












//Constructor of the cell
//We used number_Of_Microtubules = 100 in the simulation producing statistical results
//number_Of_extra_Microtubules - it is the number of microtubules that sprout asymetricaly from the MTOC just to one side
//It is intended to study the cell with assymetrical cytoskeleton
//So far, the model with assymetrical cytoskeleton was not used to get statistical results
//number_Of_extra_Microtubules = 0
Cell::Cell( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant )
{

	std::uniform_real_distribution<> distribution{ 0 , 1 };
	std::uniform_real_distribution<> distribution_angles{ 0 , 2.0 *  sim_of_Cell::PI };
	unsigned int number_of_generator = omp_get_thread_num();

	//Choosing the axis of the cell
        this->a_axis = Cell_parametres::A_AXIS;
	this->b_axis = Cell_parametres::B_AXIS;

	Vector3d center_of_nucleus( 0.0 , 0.0 , Nucleus_parametres::z_coordinate ); //1.5e-6
	this->nucleus = Nucleus( center_of_nucleus );

        //This is the control of the imput: the nucleus has to be smaller than the cell
	if( ( this->a_axis < this->nucleus.get_A_Axis() ) || ( this->b_axis < this->nucleus.get_B_Axis() ) )
	{
		cout<<"( this->a_axis < this->nucleus.get_A_Axis() ) || ( this->b_axis < this->nucleus.get_B_Axis() ) in Cell( unsigned int number_Of_Microtubules )"<<endl;
		cout<<"ERROR_ID = 976456864648645"<<endl;
		throw("");
	}
	
	//Controls, whether the same number of microtubules is used all time
	//It just checks the consistancy of the input	
	if( number_Of_Microtubules != 100 )
    	{
		cout<<"number_Of_Microtubules % 4 != 0 "<<endl;
       	 	cout<<"Cell::Cell(  unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant )"<<endl;
		cout<<"ERROR_ID = 8764356435646"<<endl;
		throw("");
    	}

    	this->number_of_microtubules = number_Of_Microtubules + number_Of_extra_Microtubules;
    	this->number_of_microtubules_extra = number_Of_extra_Microtubules;



	//////////////////////////////////////////////////////////Capture-Shrinkage////////////////////////////////////////////////////////
	//HERE EVERYTHING ABOUT IS - MICROTUBULE CATCHING WILL BE SET
	//IS_position is in x z plane - angle is determined in sim_of_Cell::IS_angle

 	double azimutal_angle = distribution_angles( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

        //The IS is constructed 
	//The IS is always placen on x-z plane 
	azimutal_angle = IS_Capture_shrinkage_param::azimutal_angle;	
	double planar_component = sin( IS_Capture_shrinkage_param::polar_angle );
	double vertical_component = cos( IS_Capture_shrinkage_param::polar_angle );
	double planar_component_x = planar_component * cos( azimutal_angle );
	double planar_component_y = planar_component * sin( azimutal_angle );


	Vector3d orientation( 0.0 , 0.0 , 0.0 );
	orientation( 0 ) = planar_component_x;
	orientation( 1 ) = planar_component_y;
	orientation( 2 ) = vertical_component;
	orientation = orientation / orientation.norm();

	//////////////////////////////////////////////////////////Capture-Shrinkage////////////////////////////////////////////////////////
	double radius_of_IS_argument = IS_Capture_shrinkage_param::radius;
	Vector3d  center_of_IS_cap_sh_fron_center = IS_Capture_shrinkage_param::z_coordinate_front * orientation;
	Vector3d  center_of_IS_cap_sh_rear_center = IS_Capture_shrinkage_param::z_coordinate_back * orientation;
	//Capture-shrinkage IS
	this->IS_Capture_Shrinkage = IS_Capt_Shrinkage( center_of_IS_cap_sh_fron_center , radius_of_IS_argument , center_of_IS_cap_sh_rear_center );




	//////////////////////////////////////////////////////////Cortical Sliding////////////////////////////////////////////////////////

    	Vector3d front_center = orientation * IS_Cortical_Sl_parameter::front_radius;
    	Vector3d rear_center = orientation * IS_Cortical_Sl_parameter::rear_radius;
    	double radius_outer = IS_Cortical_Sl_parameter::radius;
    	double radius_inner = IS_Cortical_Sl_parameter::radius_inner;
    	//Construtor of the cortical sliding IS
    	this->IS_Cortic_Sl = ISCorticalSl( front_center , rear_center , radius_outer , radius_inner );


	//Consturcting the MTOC
	//Getting the number of sprouting points
	//If the number of microtubules is bigger than MTOCparam::boundary_MTOC_points, the number of sprouting points is limited
        unsigned int number_Of_Microtubules_MTOC;
        if( number_Of_Microtubules < MTOCparam::boundary_MTOC_points )
        {
            number_Of_Microtubules_MTOC = number_Of_Microtubules;
        }
        else
        {
             number_Of_Microtubules_MTOC = MTOCparam::boundary_MTOC_points;
        }


	//The center of the MTOC is located at the x-axis
	//The distance of the MTOC from the center of coordinates
   	double z_coordinate = this->nucleus.get_B_Axis() + 3.0 / 4.0 * ( this->b_axis - this->nucleus.get_B_Axis() );
   	Vector3d center_of_MTOC( 0.0 , 0.0 , z_coordinate );
   	
   	//Constructor of the MTOC
   	this->MTOC = MTOC2( number_Of_Microtubules_MTOC , center_of_MTOC , "plane" );
   	//MTOC_simple is the object for future simulations, currently unused
   	this->mtoc = MTOC_simple( center_of_MTOC );

	//The MTOC is divided into 4 sides - each with 5 microtubule sprouting points

   	unsigned int side = 4;
   	//The number of microtubules sporuting from each side 
   	unsigned int micro_per_side = number_Of_Microtubules / side;
   	//This gives me how many microtubules is at the end of one centrioles
   	unsigned int polygon_per_side;
   	//Polygon_per_side gives me the number of polygon at every end of the centrioles
   	//It equals to the number of microtubules that will sprout from every sporuting point

   	if( micro_per_side % MTOCparam::micro_in_polygon == 0 )
   	{
       	polygon_per_side = micro_per_side / MTOCparam::micro_in_polygon;
   	}
   	else
   	{
       		polygon_per_side = micro_per_side / MTOCparam::micro_in_polygon + 1;
   	}
	//It checks the assymetrical distribution of the microtubules
   	unsigned int unfinished_polygon_number = micro_per_side % MTOCparam::micro_in_polygon;


   	Vector3d MTOC_center = this->MTOC.get_center();
   	unsigned int counter_micro = 0;
   	unsigned int polygon_counter = 0;


   	double constant = 1;
   	Vector3d center_MTOC_tmp = MTOC_center;

	//Abstract mtoc will be used in future simulations, so far unused
   	Vector3d original_orientation = this->MTOC.get_center() / this->MTOC.get_center().norm();
   	this->abstract_MTOC = ideal_MTOC( this->number_of_microtubules , original_orientation );
   	this->abstract_MTOC.set_center( this->MTOC.get_center() );
   	this->abstract_MTOC.set_orientation_center( this->MTOC.get_center() );
   	this->abstract_MTOC.set_mtoc_center( this->MTOC.get_center() );
   	this->abstract_MTOC.set_original_MTOC_center( this->MTOC.get_center() );
   	this->abstract_MTOC.set_original_MTOC_1_point( this->MTOC.get_point( 1 ) );

   	//Array holding microtubules
   	this->array_Of_Microtubules = new Microtubule[ this->number_of_microtubules ];

   	for( unsigned int strana = 0 ; strana < side ; strana ++ )
   	{
		unsigned int beads_per_side = 0;
		for( unsigned int polygon = 0 ; polygon < polygon_per_side ; polygon ++ )
		{
		    unsigned int meze;
		    if( polygon != polygon_per_side - 1 )
		    {
		        meze = MTOCparam::micro_in_polygon;
		    }
		    else
		    {
		        if( unfinished_polygon_number == 0 )
		        {
		            meze = MTOCparam::micro_in_polygon;
		        }
		        else
		        {
		            meze = unfinished_polygon_number;
		        }
		    }
		    for( unsigned int  microtubule = 0 ; microtubule < meze ; microtubule ++ )
		    {
		    	   //Two points on the MTOC are chosen and microtubule is connected to them
		           unsigned int point_number = strana * MTOCparam::micro_in_polygon  + microtubule + 1;
			   unsigned int opposite_point_number = this->MTOC.get_random_point_on_opposite_side_without_bias( point_number );
		           Vector3d MTOC_Point = this->MTOC.get_point( point_number );
		           Vector3d opposite_MTOC_Point = this->MTOC.get_point( opposite_point_number );

		           this->abstract_MTOC.set_mtoc_point(  counter_micro , MTOC_Point  );
		           Vector3d orientation = MTOC_Point - opposite_MTOC_Point;
		           orientation = orientation / orientation.norm();


		           Vector3d first_Point = opposite_MTOC_Point;
		           Vector3d second_Point = MTOC_Point;
			   std::uniform_int_distribution<int> distribution( 0 , sim_of_Cell::REDUCTION );
	    		   unsigned int number_of_generator = omp_get_thread_num();
	    		   unsigned int reduction = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

			   //reduction = 0;
			   //The number of beads of the microtubule
		           unsigned int number_of_points = sim_of_Cell::MicrotubulePoints;
		           number_of_points = sim_of_Cell::MicrotubulePoints - reduction;

		           Microtubule tmp( first_Point , second_Point , orientation , microtubule , polygon_counter , strana , point_number , opposite_point_number , number_of_points );//true

		           this->array_Of_Microtubules[ counter_micro ] = tmp;
		           this->abstract_MTOC.set_orientation( counter_micro , this->array_Of_Microtubules[ counter_micro ].getPoint( 1 ) );


		           Vector3d mtoc_point = this->MTOC.get_center() - orientation * sim_of_Cell::resting_distance;
		           this->abstract_MTOC.set_point( counter_micro , mtoc_point );
		           counter_micro = counter_micro + 1;

		    }
		    polygon_counter = polygon_counter + 1;
		}
    }
        //Abstract MTOC will be used in later  simulations, so far unused
        this->abstract_MTOC.set_true_mtoc_coordinates( this->MTOC.get_coordinates() );

}



//Returns the microtubule identified by the index
Microtubule Cell::getMicrotubule( unsigned int index )
{
	if( index >= this->number_of_microtubules )
	{
		cout<<"index >= this->number_of_microtubules in Cell::getMicrotubule( unsigned int index )"<<endl;
		throw("");
	}
	return this->array_Of_Microtubules[ index ];
}

//Returns the number of microtubules in the cell
unsigned int Cell::get_microtubule_number()
{
    return this->number_of_microtubules;
}

//Estimates the friction of the entire cytoskeleton
double Cell::get_cytoskeleton_friction()
{
     double viscosity =  sim_of_Cell::viscosity;  // 9.3 * 10.0
     double radius = 1.25e-8;
     double microtubule_lenght = sim_of_Cell::resting_distance * (double)( sim_of_Cell::MicrotubulePoints - 1 );
     double nominator = 4.0 * sim_of_Cell::PI  * microtubule_lenght; //* viscosity
     double denominator = log( microtubule_lenght / ( 2.0 * radius ) ) + 0.84;
     double effectiveFriction_micro = nominator / denominator;
     return effectiveFriction_micro * this->get_microtubule_number();

}

//Returns the capture-shrinkage
IS_Capt_Shrinkage Cell::get_IS_capture_shrinkage()
{
	IS_Capt_Shrinkage tmp;
	tmp = this->IS_Capture_Shrinkage;
	return tmp;
}



















void Cell::MidStep_3()
{

//Saves the random forces, coordinates and forces acting on microtubules
    MatrixXd* random_Forces = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* original_Coordinates = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* external_Forces = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* dynein_Forces = new MatrixXd[ this->number_of_microtubules ];

//The same for the MTOC
    MatrixXd MTOC_original_coordinates = this->MTOC.get_coordinates();
    MatrixXd MTOC_forces = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_wall_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_nucleus_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );


//Saving the coordinates of the microtubules
for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
    	original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();
}


//It has to compute forces acting on microtubule
//They differ,since the microtubule can be free, attached to cortical sliding, or capture-shrinkage

for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
        //Free microtubule
	if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() < 1 )
	{
            //external forces will contain forces from every structure in cell ( dynein , wall , MTOC )
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //Cell wall interaction
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;


            //Nucleus interaction
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC interaction
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

	    //Calulation of the MTOC force
            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );

            if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2  )
            {
                unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
                unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

                for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
                {
                    external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                    if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                    {
                          external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                    }
                }
                for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
                {
                    MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension )  ;//

                    MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
                }
            }

		external_Forces[ microtubule ] = external_Force;
	}
        else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
        {

            //external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //Cell wall interaction
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

            //Nucleus interaction
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC interaction
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );
	    //Calulation of the MTOC force
            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();


              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
            //changing force acting on the MTOC
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }

            external_Forces[ microtubule ] = external_Force;
	//Dynein forces
        MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
        this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
        external_Force = external_Force + force_dynein_IS_surface;

       external_Forces[ microtubule ] = external_Force;
    }

    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
	{
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

	    // NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );


	    if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 3 )
            {
            	this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );

            	if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2  )
            	{
            		unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            		unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

		      	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
		      	{
		          	external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
		          	if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
		          	{
		             	   external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
		          	}
		      	}
		      	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
		      	{
		         	 MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension )  ;//

		         	 MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
		      	}
           	}
            }
	    else
            {

                if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                {
            		unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            		Vector3d second_bead_force = this->MTOC_microtubule_two_points( microtubule );
              		for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              		{
                 		external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                 		MTOC_forces(3* number_MTOC_point + dimension , 0 ) = MTOC_forces( 3 * number_MTOC_point + dimension , 0 ) + (- 1) * second_bead_force( dimension );//
              		}
                }
	}
	    //Dynein
            MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
            dynein_Forces[ microtubule ] = force_dynein_IS_surface;

            external_Force = external_Force + force_dynein_IS_surface;
            external_Forces[ microtubule ] = external_Force;
    }
}

this->force_MTOC_cell_wall( MTOC_wall_force );
this->MTOC_nucleus_interaction( MTOC_nucleus_force );
MTOC_forces = MTOC_forces + MTOC_wall_force + MTOC_nucleus_force;


this->MTOC.oneStepMidStepAlgorithm_1_half( MTOC_forces );
for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ ) //this->number_of_microtubules
{
	MatrixXd randomForce = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
	this->array_Of_Microtubules[ microtubule ].oneStepMidStepAlgorithm_1_half_producing_random( randomForce , external_Forces[ microtubule ]  );
	random_Forces[ microtubule ] = randomForce;
}


MTOC_forces = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
MTOC_wall_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
MTOC_nucleus_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
this->resize_micro_with_different_first_segment_and_MTOC(  );



for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
	if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() < 1 )
	{
            //external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

	    // NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );

            if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2  )
            {
            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension )  ;//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }
            }

		external_Forces[ microtubule ] = external_Force;
	}


    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
    {

            //external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

	    // NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
            //changing force acting on the MTOC
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }

            external_Forces[ microtubule ] = external_Force;

        MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
        this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
        external_Force = external_Force + force_dynein_IS_surface;



       external_Forces[ microtubule ] = external_Force;
        //DYNEIN INTERACTION

    }

    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
	{
            //external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
            //changing force acting on the MTOC
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }


            MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );

            external_Forces[ microtubule ] = external_Force;
    }

}



//Second half of the algorithm
this->force_MTOC_cell_wall( MTOC_wall_force );
this->MTOC_nucleus_interaction( MTOC_nucleus_force );
MTOC_forces = MTOC_forces + MTOC_wall_force + MTOC_nucleus_force;
this->MTOC.oneStepMidStepAlgorithm_2_half( MTOC_forces , MTOC_original_coordinates );


for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ ) //this->number_of_microtubules
{
    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
    {
        MatrixXd force = external_Forces[ microtubule ] + dynein_Forces[ microtubule ];
        this->array_Of_Microtubules[ microtubule ].oneStepMidStepAlgorithm_2_half( original_Coordinates[ microtubule ] , random_Forces[ microtubule ] , force );
    }
    else
    {
        this->array_Of_Microtubules[ microtubule ].oneStepMidStepAlgorithm_2_half( original_Coordinates[ microtubule ] , random_Forces[ microtubule ] , external_Forces[ microtubule ] );
    }
}

    delete[] random_Forces;
    random_Forces = NULL;
    delete[] original_Coordinates;
    original_Coordinates = NULL;
    delete[] external_Forces;
    external_Forces = NULL;
    delete[] dynein_Forces;
    dynein_Forces = NULL;
}


























void Cell::Euler_algorithm()
{

//MICROTUBULES
    MatrixXd* random_Forces = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* original_Coordinates = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* external_Forces = new MatrixXd[ this->number_of_microtubules ];

//MTOC
    MatrixXd MTOC_original_coordinates = this->MTOC.get_coordinates();
    MatrixXd MTOC_forces = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_wall_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_nucleus_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );






for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
	MatrixXd randomForce = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
	this->array_Of_Microtubules[ microtubule ].getRandomForces(  sim_of_Cell::time_Step , randomForce );
	random_Forces[ microtubule ] = randomForce;

	if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() < 1 )
	{
		original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();

        //external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION

            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );


            //this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );

            if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2  )
            {
            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
            //changing force acting on the MTOC


              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension )  ;//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }
            }

		external_Forces[ microtubule ] = external_Force;
	}
    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 1 )
	{
		original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();

		//external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
		MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;


		//NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;


		//DYNEIN FORCE
            //if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 1 )
            //{
            Vector3d dyneinForceVector = this->array_Of_Microtubules[ microtubule ].dynein_force_capture_shrinkage( );
            //cout<<dyneinForceVector<<endl;
            unsigned int index_last_point = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1;
            //cout<<"index_last_point = "<<index_last_point<<endl;
            for( unsigned int i = 0 ; i < 3 ; i ++ )
            {
                external_Force( 3 * index_last_point + i , 0 ) = external_Force( 3 * index_last_point + i , 0 ) + dyneinForceVector( i );
            }



            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            //this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
            //changing force acting on the MTOC
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }
            external_Forces[ microtubule ] = external_Force;
	}
    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
    {
        original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();

        //external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            //this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
            //changing force acting on the MTOC
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }

            external_Forces[ microtubule ] = external_Force;

        MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
        //this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa( force_dynein_IS_surface );
        this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
        external_Force = external_Force + force_dynein_IS_surface;



       external_Forces[ microtubule ] = external_Force;
        //DYNEIN INTERACTION

    }

    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
	{
        original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();

        //external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            /*
            for( unsigned int i = 0 ; i < 3 ; i ++ )
            {
                wall_force_matrix( 3 * ( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 ) + i , 0 ) = 0.0;
                wall_force_matrix( 3 * ( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 ) + i , 0 ) = 0.0;
            }
            */
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            //this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
            //changing force acting on the MTOC
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }


            MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
            external_Force = external_Force + force_dynein_IS_surface;
            external_Forces[ microtubule ] = external_Force;
    }

}



//MTOC STEP MTOC STEP  MTOC STEP  MTOC STEP  MTOC STEP  MTOC STEP  MTOC STEP  MTOC STEP

this->force_MTOC_cell_wall( MTOC_wall_force );
this->MTOC_nucleus_interaction( MTOC_nucleus_force );
MatrixXd MTOC_forces_naive = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );

MTOC_forces = MTOC_forces + MTOC_wall_force + MTOC_nucleus_force;


this->MTOC.Euler_algorithm( MTOC_forces );


for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ ) //
{

	if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 2  )
	{
		continue;
	}
	else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 3 )
	{
		continue;
	}
	else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 5  )
	{
		continue;
	}
	else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 6 )
	{
		continue;
	}
	else
	{
		//if( microtubule == 2 )
		//{
			this->array_Of_Microtubules[ microtubule ].Euler_algorithm( random_Forces[ microtubule ] , external_Forces[ microtubule ] );
		//}
	}
}

    delete[] random_Forces;
    random_Forces = NULL;
    delete[] original_Coordinates;
    original_Coordinates = NULL;
    delete[] external_Forces;
    external_Forces = NULL;

}


void Cell::timeDevelopment_2( double Time_0 , double Time_1 )
{
	std::uniform_real_distribution<> distribution{0.0, 1.0};
        double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++
	//std::uniform_real_distribution<> dis(0.0, 1.0);
	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{

       		double time = step * sim_of_Cell::time_Step;
        	Cell_parametres::time = time;

        	if( step % 10 == 0)
		{
            		//cout<<"time = "<<time<<endl;
		}

        	//////////////////////////////////////

		//cout<<"number_of_generator = "<<number_of_generator<<endl;
		//unsigned int number_of_generator = omp_get_thread_num();
		//cout<<distribution(mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ])<<endl;
		this->MidStep_3();

        	//////////////////////////////////////
        	this->stepping_and_detachment_of_all_microtubule_projection_real_dynein();

        }

}



















void Cell::timeDevelopment_2( double Time_0 , double Time_1 , unsigned int number_of_Pictures )
{
	unsigned int key = 0;
        double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++

	unsigned int picture_time_step = 1;

        //cout<<"number_of_Pictures = "<<number_of_Pictures<<endl;
	if( numberOfSteps2 > number_of_Pictures )
	{
		picture_time_step = numberOfSteps2 / number_of_Pictures;
	}

	this->print_Cell_parametres();

	unsigned int printing_counter = 0;
	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{
		//cout<<"CELL step = "<<step<<endl;
       		double time = step * sim_of_Cell::time_Step;
		this->MTOC.set_time_clock( time );

        	//Cell_parametres::time = time;
        	if( numberOfSteps2 >= number_of_Pictures )
		{
			if( step % picture_time_step == 0 )
			{
				this->print_Cell( step / picture_time_step );
			}
		}



        	if( step % 1000 == 0)
		{
            		//cout<<"step % 1000 == 0"<<endl;
            		cout<<"key = "<<key<<endl;
            		cout<<"time = "<<time<<endl;
        	}


      		if( step % 50 == 0 )
		{
			//cout<<"print start"<<endl;
			this->BIG_PRINT_numerical( step , key );
                	//cout<<"print end"<<endl;
        	}


      		if( ( step %  100 ) == 0 ) //14285
		{
			if ( sim_of_Cell::density_of_dynein_motor_surface < 0.00001 )
			{	//I fear double imprecisions
				this->density_surface_dynein.change_part_of_IS_points( this->IS_Cortic_Sl );
			}
			else
			{
				this->density_surface_dynein.change_part_of_points_surface_and_IS( this->IS_Cortic_Sl );
			}
        	}



      		if( ( step %  7142 ) == 0 )
		{
			this->BIG_PRINT_numerical_micro( printing_counter , key );
			printing_counter = printing_counter + 1;
        	}

		//FUNCTIONS
        	//////////////////////////////////////
        	this->MidStep_3(); //

        	//////////////////////////////////////

        	this->stepping_and_detachment_of_all_microtubule_projection_real_dynein();


		//tohle jen kontroluje, jestli je micro chycena v IS
        	this->check_caught_micro_IS_with_real_dynein_all_micro( );
		//tohle pridava dalsi pary do mikrotubule pro IS capt shrinkage
        	this->microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro( );

		this->catch_pair_abscissa_real_dynein();


		//this->control_length_of_micro_IS();
		//this->control_length_of_micro_IS_2();
        }
}










//Saves the configuration of the cell for the visualisation
//The extended version of this function is commented in Cell_two_IScpp
void Cell::print_Cell( unsigned int index )
{

      	  char name_of_text_file [55];
	  sprintf ( name_of_text_file , "picturesVideos/textFiles/microtubule_%d.txt", index );

	  FILE *out;
	  out = fopen( name_of_text_file , "w");


	  for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
	  {

		if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 0 )
         	{

                    Vector3d tangent_last;
                    try
                    {
                        tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                    }
                    catch( int e )
                    {
                        cout<<"exception"<<endl;
                        cout<<"void void Cell::print_Cell()"<<endl;
                        throw("");
                    }

		    if( tangent_last.norm() > 2 *  sim_of_Cell::resting_distance )
		    {

		    }
		    else
		    {

	                Vector3d center = this->MTOC.get_center();
	                double parametr = 0.0;


                	for( unsigned int point = 0 ; point < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ; point ++  )
                	{
                    		parametr = parametr + 0.05;
                    		Vector3d coordinate_point = this->array_Of_Microtubules[ microtubule ].getPoint( point );
                    		fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  coordinate_point( 0 ) * 1e6 , coordinate_point( 1 ) * 1e6 , coordinate_point( 2 ) * 1e6  );
                	}
               		Vector3d tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                	Vector3d point_last = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
                	Vector3d ill_point = point_last + tangent_last / 10.0;
                	parametr = parametr + 0.05;
                	fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                	fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() );
                     	fprintf(out,"END\n" );
		   }

                }

         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
         {

		    Vector3d tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
		    if( tangent_last.norm() > 4 *  sim_of_Cell::resting_distance )
		    {
			cout<<"LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL"<<endl;
		    }
		    else
		    {
                    Vector3d center = this->MTOC.get_center();
                    //center = this->mtoc.get_position() * 1e6;
                    double parametr = 0.0;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  center( 0 )  * 1e6, center( 1 )  * 1e6, center( 2 )   * 1e6);

                    for( unsigned int point = 0 ; point < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ; point ++  )
                    {
                            parametr = parametr + 0.05;
                            Vector3d coordinate_point = this->array_Of_Microtubules[ microtubule ].getPoint( point );
                            fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  coordinate_point( 0 ) * 1e6 , coordinate_point( 1 ) * 1e6 , coordinate_point( 2 ) * 1e6  );
                    }

                    Vector3d point_last = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
                    Vector3d ill_point = point_last + tangent_last / 10.0;
                    parametr = parametr + 0.05;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                    fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() );
                    fprintf(out,"END\n" );
                    }

         }

         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
         {
                    Vector3d tangent_last;
                    try
                    {
                        tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                    }
                    catch( int e )
                    {
                        cout<<"exception"<<endl;
                        cout<<"void void Cell::print_Cell()"<<endl;
                        throw("");
                    }
		    if( tangent_last.norm() > 4 *  sim_of_Cell::resting_distance )
		    {
		    }
		    else
		    {

                    Vector3d center = this->MTOC.get_center();
                    double parametr = 0.0;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  center( 0 )  * 1e6, center( 1 )  * 1e6, center( 2 )   * 1e6);

                    for( unsigned int point = 0 ; point < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ; point ++  )
                    {
                        parametr = parametr + 0.05;
                            Vector3d coordinate_point = this->array_Of_Microtubules[ microtubule ].getPoint( point );
                        fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  coordinate_point( 0 ) * 1e6 , coordinate_point( 1 ) * 1e6 , coordinate_point( 2 ) * 1e6  );
                    }
                    Vector3d tangent_last;
                    try
                    {
                        tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                    }
                    catch( int e )
                    {
                        cout<<"exception"<<endl;
                        cout<<"void void Cell::print_Cell()"<<endl;
                        throw("");
                    }

                    Vector3d point_last = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
                    Vector3d ill_point = point_last + tangent_last / 10.0;
                    parametr = parametr + 0.05;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                    fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() );
                    fprintf(out,"END\n" );
		    }
         }
         else
	{
		cout<<endl;
		cout<<endl;
		cout<<"this->array_Of_Microtubules[ microtubule ].get_dynein_index() = "<<this->array_Of_Microtubules[ microtubule ].get_dynein_index() <<endl;
		cout<<endl;
		cout<<endl;
	}

	  }
	  fclose( out );



      //MTOC printing
	  char name_of_MTOC_file [55];
	  sprintf ( name_of_MTOC_file , "picturesVideos/textFiles/MTOC_%d.txt", index );
	  FILE *out_MTOC;
	  out_MTOC = fopen( name_of_MTOC_file , "w");
          Vector3d MTOC_point = this->MTOC.get_center();
          fprintf( out_MTOC , "%10.5f,%10.5lf,%10.5lf\n",  MTOC_point( 0 ) * 1e6 , MTOC_point( 1 ) * 1e6 , MTOC_point( 2 ) * 1e6 );
          fprintf( out_MTOC , "END\n" );


	  for( unsigned int poly_number = 1 ; poly_number <= this->MTOC.get_number_of_points() ; poly_number ++ )
	  {
                  Vector3d MTOC_point =   this->MTOC.get_point( poly_number );
		  fprintf( out_MTOC , "%10.5f,%10.5lf,%10.5lf\n",  MTOC_point( 0 ) * 1e6 , MTOC_point( 1 ) * 1e6 , MTOC_point( 2 ) * 1e6 );
		  fprintf( out_MTOC , "END\n" );
	  }
	  fclose( out_MTOC );


	  char name_of_IS_cathing_file [55];
	  sprintf ( name_of_IS_cathing_file , "picturesVideos/textFiles/IS_catching_%d.txt", index );
	  FILE *out_IS_catching;
	  out_IS_catching = fopen( name_of_IS_cathing_file , "w");

	  for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
	  {
            Vector3d IS_catching = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();
            if( IS_catching( 0 ) != 666.0 )
            {
                fprintf( out_IS_catching , "%10.5f,%10.5lf,%10.5lf\n",  IS_catching( 0 ) * 1e6 , IS_catching( 1 ) * 1e6 , IS_catching( 2 ) * 1e6 );
                fprintf( out_IS_catching , "END\n" );
            }
          }
	  fclose( out_IS_catching );










	  char IS_capture_shrinkage [55];
	  sprintf ( IS_capture_shrinkage , "picturesVideos/textFiles/IS_capture_shrinkage_%d.txt", index );
	  FILE *out_IS_capture_shrinkage;
	  out_IS_capture_shrinkage = fopen( IS_capture_shrinkage , "w");
	  Vector3d IS_capture_shrinkage_center = this->IS_Capture_Shrinkage.get_center_of_IS_front();
	  double radius_IS_Capture_Shrinkage = this->IS_Capture_Shrinkage.get_radius_of_IS();
	  Vector3d IS_capture_shrinkage_rear_point = this->IS_Capture_Shrinkage.get_center_of_IS_rear();



	  fprintf( out_IS_capture_shrinkage , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_center( 0 ) * 1e6 , IS_capture_shrinkage_center( 1 ) * 1e6 , IS_capture_shrinkage_center( 2 ) * 1e6 );
	  fprintf( out_IS_capture_shrinkage , "END\n" );
	  fprintf( out_IS_capture_shrinkage , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_rear_point( 0 ) * 1e6 , IS_capture_shrinkage_rear_point( 1 ) * 1e6 , IS_capture_shrinkage_rear_point( 2 )  * 1e6 );
	  fprintf( out_IS_capture_shrinkage , "END\n" );
	  fprintf( out_IS_capture_shrinkage , "%10.5f\n",  radius_IS_Capture_Shrinkage * 1e6 );
      fclose( out_IS_capture_shrinkage );





          char IS_cortical_sl[55];
	  sprintf ( IS_cortical_sl , "picturesVideos/textFiles/IS_cortical_sl_%d.txt", index );
	  FILE *out_IS_cortical_sl;
	  out_IS_cortical_sl = fopen( IS_cortical_sl , "w");

	  Vector3d front_center1 = this->IS_Cortic_Sl.get_center_front_of_IS();
	  Vector3d rear_center1 = this->IS_Cortic_Sl.get_center_rear_of_IS();


	  Vector3d point_on_plane = this->IS_Cortic_Sl.get_point_on_plane();
	  Vector3d axis = this->IS_Cortic_Sl.get_axis();
          Vector3d tmp_rear_point = point_on_plane + ( - 1.0 ) * axis * 2e-6;
	  double radius = IS_Cortical_Sl_parameter::radius;
          fprintf( out_IS_cortical_sl , "%10.5f,%10.5lf,%10.5lf\n",  front_center1( 0 ) * 1e6 , front_center1( 1 ) * 1e6 , front_center1( 2 ) * 1e6 );
	  fprintf( out_IS_cortical_sl , "END\n" );
	  fprintf( out_IS_cortical_sl , "%10.5f,%10.5lf,%10.5lf\n",  rear_center1( 0 ) * 1e6 , rear_center1( 1 ) * 1e6 , rear_center1( 2 )  * 1e6 );
	  fprintf( out_IS_cortical_sl , "END\n" );
	  fprintf( out_IS_cortical_sl , "%10.5f\n",  radius * 1e6 );
	  fclose( out_IS_cortical_sl );





      char MTOC_center [55];
      sprintf ( MTOC_center , "picturesVideos/textFiles/MTOC_center_%d.txt", index );
      FILE *out_MTOC_center;
      out_MTOC_center = fopen( MTOC_center , "w");
      Vector3d center_of_MTOC = this->MTOC.get_center();
      fprintf( out_MTOC_center , "%10.5f,%10.5lf,%10.5lf\n",  center_of_MTOC( 0 ) * 1e6 , center_of_MTOC( 1 ) * 1e6 , center_of_MTOC( 2 ) * 1e6 );
      fclose( out_MTOC_center );

      /////////////////////////////////////////////////////////////////////////////////////////////
           char Dynein_surface_randomly_distributed [155];
           sprintf ( Dynein_surface_randomly_distributed , "picturesVideos/textFiles/Dynein_surface_randomly_distributed_%d.txt", index );
           FILE *out_Dynein_surface_randomly_distributed;
           out_Dynein_surface_randomly_distributed = fopen( Dynein_surface_randomly_distributed , "w");
           std::map< unsigned int , std::vector<Vector3d> > mapa = this->density_surface_dynein.get_map();


           for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa.begin(); it != mapa.end(); ++it)
           {

                std::vector<Vector3d> vector_tmp = it->second;
                for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
                {
                    Vector3d vect = vector_tmp.at( i );
                    fprintf( out_Dynein_surface_randomly_distributed , "%10.5f,%10.5lf,%10.5lf\n",  vect( 0 ) * 1e6 , vect( 1 ) * 1e6 , vect( 2 ) * 1e6 );
                    fprintf( out_Dynein_surface_randomly_distributed , "END\n" );
                }
           }

           fclose( out_Dynein_surface_randomly_distributed );

          //this prints the dynein distributed in IS capture shrinkage
           char Dynein_IS_capture[155];
           sprintf ( Dynein_IS_capture , "picturesVideos/textFiles/Dynein_IS_capture_%d.txt", index );
           FILE *out_Dynein_IS_capture;
           out_Dynein_IS_capture = fopen( Dynein_IS_capture , "w");
           std::map< unsigned int , std::vector<Vector3d> > mapa_capt = this->Capture_Shrinkage_dynein.get_map();
           for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa_capt.begin(); it != mapa_capt.end(); ++it)
           {

                std::vector<Vector3d> vector_tmp = it->second;
                for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
                {
                    Vector3d vect = vector_tmp.at( i );
                    fprintf( out_Dynein_IS_capture , "%10.5f,%10.5lf,%10.5lf\n",  vect( 0 ) * 1e6 , vect( 1 ) * 1e6 , vect( 2 ) * 1e6 );
                    fprintf( out_Dynein_IS_capture , "END\n" );
                }
           }

           fclose( out_Dynein_IS_capture );



      char dynein_abscissa [155];
      sprintf ( dynein_abscissa , "picturesVideos/textFiles/dynein_abscissa_%d.txt", index );
      FILE *out_dynein_abscissa;
      out_dynein_abscissa = fopen( dynein_abscissa , "w");
      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 9 ) || ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 20 ) )
          {
                std::vector< std::pair < Vector3d ,double  > > vector_dynein_abscissa = this->array_Of_Microtubules[ micro ].get_Dynein_abscissa();
                if( vector_dynein_abscissa.size() < 1 )
                {
                    continue;
                }
                else
                {
                      for( unsigned int point = 0 ; point < vector_dynein_abscissa.size() ; point ++ )
                      {
                          std::pair < Vector3d ,double  > pair_tmp = vector_dynein_abscissa.at( point );
                          Vector3d position = std::get< 0 >( pair_tmp );
                          fprintf(out_dynein_abscissa,"%10.5f,%10.5lf,%10.5lf\n",position( 0 ) * 1e6 , position( 1 ) * 1e6 , position( 2 ) * 1e6 );
                          fprintf( out_dynein_abscissa , "END\n" );
                      }
                }
          }
      }
      fclose( out_dynein_abscissa );


      char dynein_abscissa_attachment [155];
      sprintf ( dynein_abscissa_attachment , "picturesVideos/textFiles/dynein_abscissa_attachment_%d.txt", index );
      FILE *out_dynein_abscissa_attachment;
      out_dynein_abscissa_attachment = fopen( dynein_abscissa_attachment , "w");
      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( this->array_Of_Microtubules[ micro ].get_dynein_index() == 9  )
          {
                std::vector< std::pair < Vector3d ,double  > > vector_dynein_abscissa = this->array_Of_Microtubules[ micro ].get_Dynein_abscissa();
                for( unsigned int dynein = 0 ; dynein < vector_dynein_abscissa.size() ; dynein ++ )
                {
                    std::pair < Vector3d ,double  > one_pair = vector_dynein_abscissa.at( dynein );
                    double abscissa =   std::get<1>( one_pair );
                    unsigned int lower_bead_index =  this->array_Of_Microtubules[ micro ].get_index_according_to_abscissa( abscissa );

                    Vector3d lower_point = this->array_Of_Microtubules[ micro ].getPoint( lower_bead_index );
                    Vector3d tangent = this->array_Of_Microtubules[ micro ].getTangent2( lower_bead_index );


                    double abscissa_minus_lower_bead = abscissa - ( ( double ) lower_bead_index ) * this->array_Of_Microtubules[ micro ].getRestDist();
                    Vector3d p_of_att = lower_point + ( abscissa_minus_lower_bead / tangent.norm()  ) * tangent; // this->getRestDist()
                    fprintf(out_dynein_abscissa_attachment,"%10.5f,%10.5lf,%10.5lf\n",p_of_att( 0 )*1e6,p_of_att( 1 )*1e6,p_of_att( 2 )*1e6 );
                    fprintf( out_dynein_abscissa_attachment , "END\n" );

                }
          }

      }
      fclose( out_dynein_abscissa_attachment );
	  char MTOC_SPLINES_1[55];
	  sprintf ( name_of_MTOC_file , "picturesVideos/textFiles/MTOCSpline_%d.txt", index );
	  FILE *out_MTOC_spline_1;
	  out_MTOC_spline_1 = fopen( name_of_MTOC_file , "w");


	  int number_of_segments = 10;
	  for( unsigned int poly_number = 1 ; poly_number <= this->MTOC.get_number_of_points() ; poly_number ++ )
	  {
                  Vector3d MTOC_point =   this->MTOC.get_point( poly_number );
		  Vector3d tangent_to_center = center_of_MTOC - MTOC_point;
		  Vector3d mini_tangent = tangent_to_center / (double)number_of_segments;
	          double parametr = 0.0;
		  for( int index = 0 ; index <= number_of_segments  ; index ++ )
		  {
			Vector3d tmp_point = MTOC_point + ( double ) index * mini_tangent;
			parametr = parametr + 0.05;
			fprintf(out_MTOC_spline_1,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  tmp_point( 0 ) * 1e6 , tmp_point( 1 ) * 1e6 , tmp_point( 2 ) * 1e6  );

		  }
                  fprintf(out_MTOC_spline_1,"%u\n", 1 );
                  fprintf(out_MTOC_spline_1,"END\n" );



	  }


	  for( unsigned int poly_number = 1 ; poly_number <= this->MTOC.get_number_of_points() ; poly_number ++ )
	  {
                  Vector3d MTOC_point =   this->MTOC.get_point( poly_number );
		  Vector3d next_MTOC_point( 0.0 , 0.0 , 0.0 );
		  if( poly_number < this->MTOC.get_number_of_points() )
		  {
			next_MTOC_point =  this->MTOC.get_point( poly_number + 1 );
		  }
		  else
		  {
			next_MTOC_point =  this->MTOC.get_point( 1 );
		  }

		  Vector3d tangent_to_center = next_MTOC_point - MTOC_point;
		  Vector3d mini_tangent = tangent_to_center / (double)number_of_segments;
	          double parametr = 0.0;
		  for( int index = 0 ; index <= number_of_segments  ; index ++ )
		  {
			Vector3d tmp_point = MTOC_point + ( double ) index * mini_tangent;
			parametr = parametr + 0.05;
			fprintf(out_MTOC_spline_1,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  tmp_point( 0 ) * 1e6 , tmp_point( 1 ) * 1e6 , tmp_point( 2 ) * 1e6  );

		  }
                  fprintf(out_MTOC_spline_1,"%u\n", 0 );
                  fprintf(out_MTOC_spline_1,"END\n" );



	  }
	  fclose( out_MTOC_spline_1 );




}


//Saves the basic shapes of the cell
void Cell::print_Cell_parametres( )
{
	FILE *out;
	out = fopen( "picturesVideos/cell_shapes/Cell_Parametres.txt" , "w");
	fprintf(out,"%15.12f , %15.12f\n", this->a_axis * 1e6 , this->b_axis * 1e6 );
	fclose( out );

	FILE *out2;
	out2 = fopen( "picturesVideos/cell_shapes/Nucleus_Parametres.txt" , "w");
	fprintf(out2,"%15.12f , %15.12f\n",  this->nucleus.get_A_Axis() * 1e6 , this->nucleus.get_B_Axis() * 1e6  );
	fprintf(out2,"END\n" );
	Vector3d nucleus_center = this->nucleus.get_center();
	fprintf(out2,"%15.12f , %15.12f , %15.12f", nucleus_center( 0 ) * 1e6 , nucleus_center( 1 ) * 1e6 , nucleus_center( 2 ) * 1e6 );

	fclose( out2 );
}









 //The force of the wal acting on the microtubule bead
Vector3d Cell::force_position_cell_wall_two_elipsoid( Vector3d position  )
{

    //The cell is considered to be an ellipsoid 
    //If all axii are the same - sphere
    Vector3d force( 0.0 , 0.0 , 0.0 );
    double element = 0.0;

    double BB;
    double B_multiplicator;
    if( position( 2 )  < 0 )
    {
        BB = Cell_parametres::B_AXIS_Lower + element;
        B_multiplicator = ( Cell_parametres::B_AXIS_Lower + element  ) * ( Cell_parametres::B_AXIS_Lower + element  );
    }
    else
    {
        BB =  ( Cell_parametres::B_AXIS + element );
        B_multiplicator = ( Cell_parametres::B_AXIS + element ) * ( Cell_parametres::B_AXIS + element );
    }
    double numerator = ( Cell_parametres::A_AXIS + element ) * BB;
    double A_multiplicator = ( Cell_parametres::A_AXIS + element )  * ( Cell_parametres::A_AXIS + element ) ;

    double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;
    denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
    denominator = sqrt( denominator );
    double c_multiplicator = numerator / denominator;
    if( c_multiplicator >= 1 )
    {
        return force;
    }
    else
    {
        Vector3d point_of_intersection = c_multiplicator * position;
        //distance behind the wall
        double distance_behind_wall = ( position - point_of_intersection ).norm();
        double absolut_value_of_the_force = Cell_parametres::wall_cell_k1 * ( exp( Cell_parametres::wall_cell_k2 * distance_behind_wall ) - 1.0 );
        Vector3d orientation = ( -1.0 ) * ( position - point_of_intersection ) / ( position - point_of_intersection ).norm();
        force = orientation * absolut_value_of_the_force;
        return force;
    }


}



//Simplied version of the force, the cell is taken as a sphere
Vector3d Cell::force_position_cell_wall_sphere( Vector3d position  )
{
    Vector3d force( 0.0 , 0.0 , 0.0 );
    double abs_val = position.norm();
    if( abs_val < this->a_axis )
    {

    }
    else
    {
        double distance = abs_val - this->a_axis;
        double absolut_value_of_the_force = Cell_parametres::wall_cell_k1 * ( exp( Cell_parametres::wall_cell_k2 * distance ) - 1.0 );
        Vector3d orientation = -1.0 * position / position.norm();
        force = orientation * absolut_value_of_the_force;
    }

    return force;
}










//Interaction between microtubule and the MTOC
void Cell::MTOC_microtubule_two_points_force_with_bending( Vector3d& first_point_force , Vector3d& second_point_force , Vector3d& bending_opposite_MTOC_point , Vector3d& bending_MTOC_point , unsigned int microtubule_number )
{
	//https://www.rieger.uni-saarland.de/Paper/supplementary_final_hornak_rieger_2020.pdf
    if( microtubule_number >= this->number_of_microtubules )
    {
        cout<<"MTOC_microtubule_two_points_force_with_bending( MatrixXd force , unsigned int microtubule_number )"<<endl;
        cout<<"microtubule_number >= this->number_of_microtubules"<<endl;
        cout<<"microtubule_number = "<<microtubule_number<<endl;
        cout<<"ERROR_ID = 666719816168645"<<endl;
        throw("");
    }

    if( this->array_Of_Microtubules[ microtubule_number ].getNumberOfPoints() >= 2 )
    {

        unsigned int mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_point();
        unsigned int opposite_mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_opposite_point();
        Vector3d MTOC_opposite_point = this->MTOC.get_point( opposite_mtoc_point_index );
        Vector3d first_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint(0);

        double distance_1 = ( MTOC_opposite_point - first_point_micro ).norm();
        if( distance_1 == 0 )
        {

        }
        else
        {
            Vector3d orientation_1 = ( MTOC_opposite_point - first_point_micro );
            orientation_1 = orientation_1 / orientation_1.norm();
            double force_1_abs = distance_1 * MTOCparam::MTOC2_point_micro_kappa;
            Vector3d force_1 = force_1_abs * orientation_1;
            first_point_force = force_1;
        }



        Vector3d MTOC_point = this->MTOC.get_point( mtoc_point_index );
        Vector3d second_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint( 1 );

        double distance_2 = ( MTOC_point - second_point_micro ).norm();
        if( distance_2 == 0 )
        {

        }
        else
        {
            Vector3d orientation_2 = ( MTOC_point - second_point_micro );
            orientation_2 = orientation_2 / orientation_2.norm();
            double force_2_abs = distance_2 * MTOCparam::MTOC2_point_micro_kappa;
            Vector3d force_2 = force_2_abs * orientation_2;
            second_point_force = force_2;
        }

	double stiffness_coeffcient = 1.0;

        //add bending force
        Vector3d bending_force_one_bead_micro( 0.0 , 0.0 , 0.0 );

        Vector3d mic_first_segment = this->array_Of_Microtubules[ microtubule_number ].getTangent2( 0 );
        Vector3d MTOC_segment = MTOC_point - MTOC_opposite_point;

        Vector3d first = mic_first_segment;
	Vector3d second = MTOC_segment;


        double magThisVec = first.norm();
	first = first * ( 1.0 / magThisVec );

	double magNextVec = second.norm();
	second = second * ( 1.0 / magNextVec );		//ATTENTION - here, it will be divided by magNextVec - the last tangent

	bending_force_one_bead_micro( 0 ) = ( first( 0 ) - second( 0 ) * first.dot( second ) ) / magNextVec;
	bending_force_one_bead_micro( 1 ) = ( first( 1 ) - second( 1 ) * first.dot( second ) ) / magNextVec;
	bending_force_one_bead_micro( 2 ) = ( first( 2 ) - second( 2 ) * first.dot( second ) ) / magNextVec;
        bending_force_one_bead_micro = ( -1.0 ) * sim_of_Cell::k_bending_analytical / magThisVec * bending_force_one_bead_micro * stiffness_coeffcient;
        second_point_force = second_point_force + bending_force_one_bead_micro;


        Vector3d bending_force_MTOC_point( 0.0 , 0.0 , 0.0 );
        bending_force_MTOC_point( 0 ) = ( - second( 0 ) + first( 0 ) * first.dot( second ) ) / magThisVec;
        bending_force_MTOC_point( 1 ) = ( - second( 1 ) + first( 1 ) * first.dot( second ) ) / magThisVec;
	bending_force_MTOC_point( 2 ) = ( - second( 2 ) + first( 2 ) * first.dot( second ) ) / magThisVec;
        bending_force_MTOC_point = ( -1.0 ) * sim_of_Cell::k_bending_analytical / magNextVec * bending_force_MTOC_point * stiffness_coeffcient;
        bending_MTOC_point = bending_force_MTOC_point;

        Vector3d central_point_force = ( -1.0 ) * bending_force_one_bead_micro + ( -1.0 ) * bending_force_MTOC_point;
        first_point_force = first_point_force + central_point_force;

        bending_opposite_MTOC_point = central_point_force;


    }
    else if( this->array_Of_Microtubules[ microtubule_number ].getNumberOfPoints() == 1 )
    {
        unsigned int opposite_mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_opposite_point();
        Vector3d MTOC_opposite_point = this->MTOC.get_point( opposite_mtoc_point_index );
        Vector3d first_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint(0);

        double distance_1 = ( MTOC_opposite_point - first_point_micro ).norm();
        if( distance_1 == 0 )
        {

        }
        else
        {
            Vector3d orientation_1 = ( MTOC_opposite_point - first_point_micro );
            orientation_1 = orientation_1 / orientation_1.norm();
            double force_1_abs = distance_1 * MTOCparam::MTOC2_point_micro_kappa;
            Vector3d force_1 = force_1_abs * orientation_1;
            first_point_force = force_1;

        }
    }




}

//Simplified calculation of the force between the microtubule and MTOC
Vector3d Cell::MTOC_microtubule_two_points(  unsigned int microtubule_number )
{
    if( microtubule_number >= this->number_of_microtubules )
    {
        cout<<"MTOC_microtubule_two_points( MatrixXd force , unsigned int microtubule_number )"<<endl;
        cout<<"microtubule_number >= this->number_of_microtubules"<<endl;
        cout<<"microtubule_number = "<<microtubule_number<<endl;
        cout<<"ERROR_ID = 6667198641546865645"<<endl;
        throw("");
    }



    unsigned int mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_point();
    Vector3d MTOC_point = this->MTOC.get_point( mtoc_point_index );
    Vector3d second_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint( 1 );

    double distance_2 = ( MTOC_point - second_point_micro ).norm();
    if( distance_2 == 0 )
    {
	Vector3d force( 0.0 , 0.0 , 0.0 );
	return force;
    }
    else
    {
        Vector3d orientation_2 = ( MTOC_point - second_point_micro );
        orientation_2 = orientation_2 / orientation_2.norm();
        double force_2_abs = distance_2 * MTOCparam::MTOC2_point_micro_kappa;
        Vector3d force_2 = force_2_abs * orientation_2;
        return force_2;
    }



}




//Computes the force of the nucleus acting on the entire microtubule given by index
void Cell::force_whole_microtubule_nucleus( MatrixXd &force , unsigned int index )
{
	if( index >= this->number_of_microtubules )
	{
		cout<<"index >= this->number_of_microtubules in Cell::force_whole_microtubule_nucleus( MatrixXd &force , unsigned int index )"<<endl;
		cout<<"index = "<<index<<endl;
		throw("");
	}
	Microtubule temporary_microtubule = this->getMicrotubule( index );

	if( force.rows() != 3 * temporary_microtubule.getNumberOfPoints()  )
	{
		cout<<"force.rows() != 3 * temporary_microtubule.getNumberOfPoints() in force_whole_microtubule_nucleus"<<endl;
		cout<<"index = "<<index<<endl;
		throw("");
	}

	for( unsigned int bead_index = 0 ; bead_index < temporary_microtubule.getNumberOfPoints() ; bead_index ++ )
	{
		Vector3d position = temporary_microtubule.getPoint( bead_index );
		Vector3d force_one_bead = this->nucleus.force_position_nucleus_wall( position );

		for( unsigned int index_dimension = 0 ; index_dimension < 3 ; index_dimension ++ )
		{
			force( 3 * bead_index + index_dimension , 0 ) = force_one_bead( index_dimension );
		}
	}

}


//Computes the force acting on one point of the MTOC - uses force_position_cell_wall_two_elipsoid
Vector3d Cell::MTOC_point_wall_interaction( Vector3d position )
{
     return force_position_cell_wall_two_elipsoid( position );
}



//Computes the force of the wall on  the MTOC
void Cell::force_MTOC_cell_wall( MatrixXd &force_MTOC )
{

    if( force_MTOC.rows() != 3 * ( this->MTOC.get_number_of_points() + 1 ) )
    {
        cout<<"force_MTOC.rows() != this->MTOC.get_number_of_points() * 3"<<endl;
        cout<<"Cell::force_MTOC_cell_wall( MatrixXd &force_MTOC )"<<endl;
        cout<<"ERROR_ID Cell 4646116911546945"<<endl;
        throw("");
    }
    //Center	
    Vector3d center = this->MTOC.get_center();
    Vector3d force = this->MTOC_point_wall_interaction( center );

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        force_MTOC( dimension , 0 ) = force( dimension );
    }
}


//Force of the MTOC acting on all MTOC points
void Cell::MTOC_nucleus_interaction( MatrixXd &force_on_MTOC )
{

    if( ( force_on_MTOC.rows() != ( this->MTOC.get_number_of_points() + 1 ) * 3 ) || ( force_on_MTOC.cols() != 1 ) )
    {
        cout<<"( force_on_MTOC.rows() != ( this->MTOC.get_number_of_points() + 1 ) * 3 ) || ( force_on_MTOC.cols() != 1 )"<<endl;
        cout<<" Cell::MTOC_nucleus_interaction( MatrixXd &force_on_MTOC ) "<<endl;
        cout<<"CELL ERROR_ID = 68464186441684"<<endl;
        throw("");
    }

    Vector3d center_position = this->MTOC.get_center();
    Vector3d force_center = this->nucleus.force_position_nucleus_wall( center_position );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        force_on_MTOC( dimension , 0 ) = force_center( dimension , 0 );
    }

    for( unsigned int bod = 1 ; bod <= this->MTOC.get_number_of_points() ; bod ++ )
    {
        Vector3d position_point = this->MTOC.get_point( bod );
        Vector3d force_on_point = this->nucleus.force_position_nucleus_wall( position_point );
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            force_on_MTOC( 3 * bod + dimension , 0 ) = force_on_point( dimension , 0 );
        }
    }


}








//Sets the map holding the cortical sliding dyneins
void Cell::set_density_surface_dynein( Surface& density )
{
    this->density_surface_dynein = density;
}






//This function resizes the cytoskeleton after every step
//The projection operator project forces on hyperplane so that they are perpendicular to the constraints
//Still, the lenghts of the microtubules can change due to numerical imprecision
//This function corrects the lengths
void Cell::stepping_and_detachment_of_all_microtubule_projection_real_dynein()
{

    //resizing of the MTOC
    this->MTOC.resize_from_originals(  );

//MICROTUBULES

    for( unsigned int micro_index = 0 ; micro_index < this->number_of_microtubules ; micro_index ++ )
    {
        if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 20  )
        {
            //Detached dyneins
	    std::vector<Vector3d> motors_for_surface = this->array_Of_Microtubules[ micro_index ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2();
	    IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
	    //Creation of the new point
            std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( motors_for_surface.size() , tmp );
            this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
        }
        else if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 9  )
        {
            //Resizing of the micro
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
            //Stepping of dynein
            std::vector<Vector3d> motors_for_surface = this->array_Of_Microtubules[ micro_index ].stepping_detach_real_dynein_abscissa_projection();



	    unsigned int number_of_motors_for_IS_1 = 0;
            unsigned int number_in = 0;
            unsigned int number_out = 0;
            //Controls whether dyneins belong to the cortical sliding IS, or to the membrane
    	    std::vector<Vector3d> vector_to_project = this->IS_Cortic_Sl.control_IS_Cortical_Sliding_points2( motors_for_surface , number_in , number_out );


	    //Projection and adding of point on the surface	
 	    this->project_and_add_points_to_surface( vector_to_project );

	    if( number_in > 0 )
	    {
		std::vector<Vector3d> motors_for_surface = this->IS_Cortic_Sl.create_IS_Cortical_Sliding_points( number_in );
                this->project_and_add_points_to_surface( motors_for_surface );
	    }

	    if( number_out > 0 )
	    {
		std::vector<Vector3d> positions_outside_IS = this->density_surface_dynein.create_Cortical_Sliding_points_outside_IS( number_out , this->IS_Cortic_Sl );
                this->project_and_add_points_to_surface( positions_outside_IS );
	    }

        }
        else
        {
            //resizing of a free microtubule
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
        }
    }





}


















//This is the function for control, whether the microtubule attaches to capture-shrinkage dynein
void Cell::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )
{
    //This is the function for control, whether the microtubule 0 or 9 is caught
    //The microtubule has to gain the point, where it is depolimerized
    //It is the original anchor point of dynein 


    std::uniform_real_distribution<> distribution{ 0 , 1 };

    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
    {
	cout<<"void Cell::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )"<<endl;
	cout<<" this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 "<<endl;
	unsigned int ERROR_ID = 978161641;
	cout<<"ERROR_ID = "<<ERROR_ID<<endl;
	throw("");
    }

    unsigned int dynein_index_to_remember;
    //Only the last segment can attach
    unsigned int lower_index_tmp = 1;
    if( dynamic_instability::capture_shrinkage_switch == true )
    {
	lower_index_tmp = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2;
    }	 	

    for( unsigned int segment_id = lower_index_tmp ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 ; segment_id ++ )
    {
	if( ( ( segment_id == 1 ) || ( segment_id == 2 ) ) && ( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 5 )   )
	{
		continue;
	}
	//Bead position and tangent
        Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d tangent;
        try
        {
            tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );
        }
        catch( int e )
        {
            cout<<"exception"<<endl;
            cout<<"void Cell::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )"<<endl;
            throw("");
        }
        
        //Checks whether segment passes nearby capture-shrinkage IS
        bool answer = this->IS_Capture_Shrinkage.check_IS_capture_segment_control( bead_position , tangent );
        
        //Checks the attachment of dyneins
        if( answer == true  )
        {
            //Get the dyneins in the IS
            std::vector<Vector3d> dynein_points = this->Capture_Shrinkage_dynein.get_dynein_points( 1 );
            std::vector< Vector3d > replace_motors;


	    //Checks whether dynein attaches
            for( unsigned int dynein_index = 0 ; dynein_index < dynein_points.size() ; dynein_index ++ )
            {
                Vector3d motor_position = dynein_points[ dynein_index ];
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                //Computes the closest point and the distance
                double distance = distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );
		unsigned int number_of_generator = omp_get_thread_num();
    		double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

                if( probability < Dynein_real::calculate_attachment_probability_per_time_step( distance ) )
                {
			dynein_index_to_remember = dynein_index;
			//erasing all dynein from nine
                        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9  )
                        {
                            std::vector< Vector3d  > points_positions = this->array_Of_Microtubules[ microtubule ].get_Dynein_points_and_erase();
                            this->project_and_add_points_to_surface( points_positions );
                        }
                        
			Vector3d position = motor_position;/////////////////////////////////////////////////////////
                        Vector3d caught_position = this->project_point_on_surface( position );

			//This just controls whether the microtubule is not prolonged more than a treshold
			//The microtubule can prolong, since the + end is placed on the surface of the cell
			if( segment_id == this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2  )
			{
				Vector3d first_point = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
				double distance_new_tangent = ( first_point - caught_position ).norm();
				double treshold_tmp = 1e-6;
				double treshold_tmp_2 = this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule_outside_MTOC() / 5.0;
				if( treshold_tmp_2 < treshold_tmp ) 
				{
					treshold_tmp = treshold_tmp_2;
				}
				if(  distance_new_tangent > this->array_Of_Microtubules[ microtubule ].get_last_Tangent( ).norm() + treshold_tmp  )
				{
					continue;
				}
				else
				{

				}
			}

			//Adds point where the microtubule depolimerizes
                        this->array_Of_Microtubules[ microtubule ].set_IS_position_catching( caught_position  );
                        unsigned int number_of_beads = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1;

                        //substraction of the beads

                        for( unsigned int sub_id = number_of_beads ; sub_id > segment_id ; sub_id -- )
                        {
                        	this->array_Of_Microtubules[ microtubule ].Substract_one_Bead();
                        }
                        this->array_Of_Microtubules[ microtubule ].add_one_final_Bead_catching_position();

			this->array_Of_Microtubules[ microtubule ].set_dynein_index( 20 );
			this->array_Of_Microtubules[ microtubule ].set_growth_index( 0 );
                        this->array_Of_Microtubules[ microtubule ].set_lenght_after_catching( this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule() );
			this->array_Of_Microtubules[ microtubule ].set_lenght_of_tangents("void Cell::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )");
			break;
		}
	    }


	    //break;
	}
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20  )
    {
    	this->array_Of_Microtubules[ microtubule ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple_after_catching();
    }


    //Attachment of the dynein point
    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20  )
    {

    	std::vector<Vector3d> dynein_points = this->Capture_Shrinkage_dynein.get_dynein_points( 1 );
    	
        Vector3d motor_position = dynein_points[ dynein_index_to_remember ];


    	unsigned int bead_index = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2;
    	Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( bead_index );
    	Vector3d tangent;
    	try
    	{
        	tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( bead_index );
    	}
    	catch( int e )
    	{
        	cout<<"exception"<<endl;
        	cout<<"void Cell::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )"<<endl;
        	throw("");
    	}

	//The dynein is attached and anchor point and the abscissa are created
   	Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
   	distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );
   	double distance_lower_bead = 0;
   	for( unsigned int i = 0 ; i < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 ; i ++ )
   	{
   		distance_lower_bead = distance_lower_bead + this->array_Of_Microtubules[ microtubule ].getTangent2( i ).norm();
   	}

   	double distance_point_attachment = ( bead_position - cl_p_of_seg ).norm();
   	double abscissa = ( distance_lower_bead + distance_point_attachment );
   	if( distance_point_attachment == tangent.norm() )
   	{
   		abscissa = abscissa - 1e-8;
   	}

   	if( abscissa > this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule() )
   	{
   		throw("");
  	 }

  	std::pair < Vector3d , double > tmp_pair( cl_p_of_seg , abscissa );
  	//Pair is added
   	this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );


    	dynein_points.erase (dynein_points.begin() + dynein_index_to_remember );
    	//The rest of the point is returned to IS
    	this->Capture_Shrinkage_dynein.set_dynein_points( 1 , dynein_points );

    }


}




//Checks whether a new microtubule is attached to capture-shrinkage dynein
void Cell::check_caught_micro_IS_with_real_dynein_all_micro( )
{
    Vector3d center_of_MTOC = this->MTOC.get_center();
    Vector3d center_of_IS = this->IS_Capture_Shrinkage.get_center_of_IS_front();
    double distance = ( center_of_MTOC - center_of_IS ).norm();

    if( distance > 0.0e-6 )
    {
        for( unsigned int micro_id = 0 ; micro_id < this->get_microtubule_number() ; micro_id ++ )
        {
            if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 20 )
            {
                this->check_caught_micro_IS_with_real_dynein( micro_id );
            }
        }
    }
}



//Adds new dynein to the microtubule already attached to the first IS
//The function is the same as void Cell_two_IS::microtubule_catch_pair_abscissa_real_dynein_in_IS_2( unsigned int microtubule )
//Comments can be found in mentioned function
//The only difference is that it dels with microtubules attached to the first IS
void Cell::microtubule_catch_pair_abscissa_real_dynein_in_IS( unsigned int microtubule )
{

    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();
    unsigned int lower_index_MAIN = 1;
    if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 6 )
    {
	lower_index_MAIN = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 4;
    }	

    for( unsigned int segment_id = lower_index_MAIN ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 ; segment_id ++ )
    {


        Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d tangent;
        try
        {
            tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );
        }
        catch( int e )
        {
            cout<<"exception"<<endl;
            cout<<"void void Microtubule::set_lenght_of_tangents()"<<endl;
            throw("");
        }


        bool answer = this->IS_Capture_Shrinkage.check_IS_capture_segment_control( bead_position , tangent );
        if( ( answer == true ) && ( this->Capture_Shrinkage_dynein.get_dynein_points( 1 ).size() > 0 ) )
        {

            std::vector<Vector3d> replace_motors;
            std::vector<Vector3d> motors_to_attach;
            std::vector<Vector3d> dynein_points = this->Capture_Shrinkage_dynein.get_dynein_points( 1 );

            for( unsigned int dynein_index = 0 ; dynein_index < dynein_points.size() ; dynein_index ++ )
            {

                Vector3d motor_position = dynein_points[ dynein_index ];
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );

		double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

                if( probability < Dynein_real::calculate_attachment_probability_per_time_step( distance ) )
                {
                    motors_to_attach.push_back( motor_position );
                }
                else
                {
                    replace_motors.push_back( motor_position );
                }
            }


            for( unsigned int counter_id = 0 ; counter_id < motors_to_attach.size() ; counter_id ++ )
            {
                Vector3d motor_position = motors_to_attach[ counter_id ];
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );

                double distance_lower_bead = this->array_Of_Microtubules[ microtubule ].get_distance_to_lower_bead_with_index( segment_id );
                double distance_point_attachment = ( bead_position - cl_p_of_seg ).norm();
                double abscissa = ( distance_lower_bead + distance_point_attachment );

                if( abs( abscissa - this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule() ) < Dynein_real::step / 2.0 )
                {
                    abscissa = abscissa - Dynein_real::step;
                }

                std::pair < Vector3d , double > tmp_pair( cl_p_of_seg , abscissa );
                this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
            }


            if( replace_motors.size() == 0 )
            {
                this->Capture_Shrinkage_dynein.erase_vector_points_with_key( 1 );
            }
	    else
	    {

                this->Capture_Shrinkage_dynein.set_dynein_points( 1 , replace_motors );
            }
        }
    }
}

//Adds new capture-shrinkage dyneins
void Cell::microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro( )
{
    for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++ )
    {
            if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
            {
                this->microtubule_catch_pair_abscissa_real_dynein_in_IS( microtubule );
            }
    }

}

//Checks the length of microtubule because of the numerical inprecisions during attachment
void Cell::control_length_of_micro_IS()
{
    //It detaches the microtubules


    for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++ )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
        {

            double original_lenght = this->array_Of_Microtubules[ microtubule ].get_lenght_after_catching();
            double current_lenght = this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule();

            if( current_lenght < original_lenght )
            {
                this->array_Of_Microtubules[ microtubule ].set_lenght_after_catching( current_lenght );
            }

            if( current_lenght > original_lenght * IS_Capture_shrinkage_param::procentage_constant )
            {
		//cout<<"Cell::control current_lenght > original_lenght * IS_Capture_shrinkage_param::procentage_constant BBBBBBBBBBBBBBB"<<endl;

                std::vector<Vector3d> erased_vectors = this->array_Of_Microtubules[ microtubule ].get_dynein_points_in_IS_and_erase();
		IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
                std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( erased_vectors.size() , tmp );
                this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
                this->array_Of_Microtubules[ microtubule ].set_dynein_index( 0 );

		Vector3d posledni_tangent = this->array_Of_Microtubules[ microtubule ].get_last_Tangent();
		Vector3d position = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
		position = position + posledni_tangent / posledni_tangent.norm() * sim_of_Cell::resting_distance;
		this->array_Of_Microtubules[ microtubule ].setPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 , position );
		//cout<<"this->array_Of_Microtubules[ microtubule ].get_last_Tangent().norm() = "<<this->array_Of_Microtubules[ microtubule ].get_last_Tangent().norm()<<endl;
		this->array_Of_Microtubules[ microtubule ].set_lenght_of_tangents("void Cell::control_length_of_micro_IS()");


            }

            if( this->array_Of_Microtubules[ microtubule ].get_last_Tangent().norm() > 4 * sim_of_Cell::resting_distance )
            {
		//cout<<"Cell::control this->array_Of_Microtubules[ microtubule ].get_last_Tangent().norm() > 4 * sim_of_Cell::resting_distance"<<endl;
                std::vector<Vector3d> erased_vectors = this->array_Of_Microtubules[ microtubule ].get_dynein_points_in_IS_and_erase();
		IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
                std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( erased_vectors.size() , tmp );
                this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
                this->array_Of_Microtubules[ microtubule ].set_dynein_index( 0 );

		Vector3d posledni_tangent = this->array_Of_Microtubules[ microtubule ].get_last_Tangent();
		Vector3d position = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
		position = position + posledni_tangent / posledni_tangent.norm() * sim_of_Cell::resting_distance;
		this->array_Of_Microtubules[ microtubule ].setPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 , position );
		//cout<<"this->array_Of_Microtubules[ microtubule ].get_last_Tangent().norm() = "<<this->array_Of_Microtubules[ microtubule ].get_last_Tangent().norm()<<endl;
		this->array_Of_Microtubules[ microtubule ].set_lenght_of_tangents("void Cell::control_length_of_micro_IS()");

            }





        }
    }

}





void Cell::control_length_of_micro_IS_2()
{
	//tohle se tyka jen microtubule 20 - ohyb, odtrhnuti


    for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++ )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
        {
            std::vector<Vector3d> erased_vectors = this->array_Of_Microtubules[ microtubule ].control_length_of_micro_IS_2();
	    if( erased_vectors.size() >= 1 )
	    {
	    	IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
            	std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( erased_vectors.size() , tmp );
            	this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
            }

        }
    }

}







//It takes the point and projects it on plasma membrane
//Plasma membrane is considered as an ellipse - sphere in current simulations
Vector3d Cell::project_point_on_surface( Vector3d position )
{
    double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;
    double B_multiplicator;
    double numerator;
    if( position( 2 ) > 0 )
    {
        B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        numerator = Cell_parametres::A_AXIS * Cell_parametres::B_AXIS;
    }
    else
    {
        B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        numerator = Cell_parametres::A_AXIS * Cell_parametres::B_AXIS_Lower;
    }

    double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;
    denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
    denominator = sqrt( denominator );
    double c_multiplicator = numerator / denominator;

    Vector3d new_position = position * c_multiplicator;
    return new_position;
}



//It returns the id of the compartment containing the point
unsigned int Cell::get_dynein_compartment_id( Vector3d position )
{
    Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );
    Vector3d position_to_get_index = position + add_to_get_segment;

    unsigned int x_index = position_to_get_index( 0 ) / this->density_surface_dynein.get_X_width();
    unsigned int y_index = position_to_get_index( 1 ) / this->density_surface_dynein.get_Y_width();
    unsigned int z_index = position_to_get_index( 2 ) / this->density_surface_dynein.get_Z_width();

    unsigned int neccessary_dimension = this->density_surface_dynein.get_neccessary_dimension();
    unsigned int ID_map_index = x_index * neccessary_dimension * neccessary_dimension + y_index * neccessary_dimension + z_index;
    return ID_map_index;
}




//Projects the points on plasma membrane and adds them to the map
void Cell::project_and_add_points_to_surface( std::vector< Vector3d > points )
{
    for( unsigned int point_id = 0 ; point_id < points.size() ; point_id ++ )
    {
        Vector3d point_position = points.at( point_id );
        Vector3d point_position_new = this->project_point_on_surface( point_position );
        unsigned int map_key = get_dynein_compartment_id( point_position_new );
        this->density_surface_dynein.add_dynein_point( map_key , point_position_new );
    }

}

//Adds points to the map
void  Cell::add_points_to_surface( std::vector< Vector3d > points )
{
    for( unsigned int point_id = 0 ; point_id < points.size() ; point_id ++ )
    {
        Vector3d point_position = points.at( point_id );
        unsigned int map_key = get_dynein_compartment_id( point_position );
        this->density_surface_dynein.add_dynein_point( map_key , point_position );
    }
}




//Resizing of the cytoskeleon due to numerical imprecisions
void Cell::resize_micro_with_different_first_segment_and_MTOC(  )
{
    this->MTOC.resize_from_originals();

    for( unsigned int micro_index = 0 ; micro_index < this->number_of_microtubules ; micro_index ++ )
    {
        if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 20  )
        {
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
        }
        else
        {
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
        }
    }



}




//Cortical sliding dynein attaches on microtubule
void Cell::catch_pair_abscissa_real_dynein()
{
    for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 0 )
        {
            this->microtubule_catch_pair_abscissa_real_dynein( microtubule );
        }
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
        {
            this->microtubule_catch_pair_abscissa_real_dynein( microtubule );
        }

    }
}





//Returns the all dyneins in the compartment determined by ID_map_index
std::vector<Vector3d> Cell::get_dynein_in_compartment( unsigned int ID_map_index )
{

    std::vector<Vector3d> returned_points;
    std::vector<Vector3d> dynein_on_surface_points = this->density_surface_dynein.get_dynein_points( ID_map_index );
    returned_points = dynein_on_surface_points;
    return returned_points;
}












//Get number of cortical sliding dyneins attached to microtubules
unsigned int Cell::get_number_of_dynein_motors_cortical_sliding_micro()
{

    unsigned int sum_of_motor_in_micro = 0;
    for( unsigned int motor_id = 0 ; motor_id < this->number_of_microtubules; motor_id ++ )
    {
        if( this->array_Of_Microtubules[ motor_id ].get_dynein_index() == 9 )
        {
            unsigned int number_of_dynein_in_micro = this->array_Of_Microtubules[ motor_id ].get_Dynein_abscissa().size();
            sum_of_motor_in_micro = sum_of_motor_in_micro + number_of_dynein_in_micro;
        }
    }

    return sum_of_motor_in_micro;
}



























//Gets the number of attached capture-shrinkage dyneins
unsigned int Cell::get_number_of_dynein_motors_in_micro_capture()
{
    unsigned int number_of_motors = this->Capture_Shrinkage_dynein.get_dynein_motors_number();

    unsigned int sum_of_motor_in_micro = 0;
    unsigned int number_of_20_micro = 0;
    //Goes throught the microtubules and adds capture-shrinkage dyneins
    for( unsigned int motor_id = 0 ; motor_id < this->number_of_microtubules; motor_id ++ )
    {
        if( this->array_Of_Microtubules[ motor_id ].get_dynein_index() == 20 )
        {
            unsigned int number_of_dynein_in_micro = this->array_Of_Microtubules[ motor_id ].get_Dynein_abscissa().size();
            sum_of_motor_in_micro = sum_of_motor_in_micro + number_of_dynein_in_micro;
	    number_of_20_micro = number_of_20_micro + 1;
        }
    }
    return sum_of_motor_in_micro;
}





//Adds new dynein to microtubule
void Cell::microtubule_catch_pair_abscissa_real_dynein( unsigned int microtubule )
{

    std::uniform_real_distribution<> distribution{ 0 , 1 };
    double lenght_micro = this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule();
    //checks whether segment pass throught the compartment containing the dynein
    for( unsigned int segment_id = 1 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  - 1 ; segment_id ++ )
    {

        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d lower_bead_position_2 = this->project_point_on_surface( lower_bead_position );
        Vector3d upper_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id + 1 );
        Vector3d upper_bead_position_2 = this->project_point_on_surface( upper_bead_position );
        Vector3d tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );
        Vector3d tangent_2 = upper_bead_position_2 - lower_bead_position_2;
        Vector3d mini_segment = tangent_2 / ( double ) Dynein::number_of_segment_steps;

	//Here are saved the ids of all compartments
        std::vector< unsigned int > compartment_IDs;
        unsigned int last_index = 0;
        for( unsigned int counter = 0 ; counter < Dynein::number_of_segment_steps ; counter ++ )
        {
            Vector3d tmp_position = lower_bead_position_2 + ( double ) counter * mini_segment;
            unsigned int compartment_id = this->density_surface_dynein.get_dynein_compartment_id_projected( tmp_position );
            if( counter == 0 )
            {

                if( compartment_id != 0 )
                {
                    compartment_IDs.push_back( compartment_id );
                }
            }
            else
            {
                if( compartment_id != 0 )
                {
                    if( compartment_id != last_index )
                    {
                        compartment_IDs.push_back( compartment_id );
                    }
                }
            }
            last_index = compartment_id;
        }


        for( unsigned int comp_index = 0 ; comp_index < compartment_IDs.size() ; comp_index ++ )
        {
            unsigned int compartment_id = compartment_IDs.at( comp_index );
	    //Returns dyneins in the compartment with the ID
            std::vector<Vector3d> dynein_points = this->get_dynein_in_compartment( compartment_id );
            //These points will be returned back
            std::vector< Vector3d > replace_motors;
            for( unsigned int bod_id  = 0 ; bod_id < dynein_points.size() ; bod_id ++ )
            {
                Vector3d motor_position = dynein_points.at( bod_id );
                //Computation of the closest point and the distance
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                //Distance
                double distance = distance_point_segment( lower_bead_position , tangent , motor_position , cl_p_of_seg );
                unsigned int number_of_generator = omp_get_thread_num();
		double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
		//The attachment of the point
                    if( probability < Dynein_real::calculate_attachment_probability_per_time_step( distance ) )
                    {
                        double distance_lower_bead = 0;
                        for( unsigned int i = 0 ; i < segment_id ; i ++ )
                        {
                            distance_lower_bead = distance_lower_bead + this->array_Of_Microtubules[ microtubule ].getTangent2( i ).norm();
                        }

                        double distance_point_attachment = ( lower_bead_position - cl_p_of_seg ).norm();
                        //Computation of the abscissa
                        double abscissa = ( distance_lower_bead + distance_point_attachment ) ;

                        if( abscissa >= lenght_micro - 8.0e-9 )
                        {
                            replace_motors.push_back( motor_position );
                            continue;
                        }
                        else
                        {
			    //Setting the position of the anchor point
                            Vector3d new_motor_position = cl_p_of_seg;
                            //Creating the pair
                            std::pair < Vector3d , double > tmp_pair( new_motor_position , abscissa );
                            //Adding the dynein to microtubule
                            this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
                            this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );
                        }

                    }
                    else
                    {	
                    	//If the dynein does not attach
                        replace_motors.push_back( motor_position );
                    }

            }
            //Points are returned
            this->density_surface_dynein.set_dynein_points( compartment_id , replace_motors );
        }



    }

}











//Saves the time and the center of the MTOC
void Cell::BIG_PRINT( unsigned int time_Step )
{
    this->print_time( time_Step );
    this->print_center_MTOC( );
}




//Saves the center of the MTOC, number of attached capture-shrinkage and cortical sliding dyneins into files with the name containing the id of the process
void Cell::BIG_PRINT_numerical( unsigned int time_Step , unsigned int key )
{
    //MTOC center
    this->print_center_MTOC_numerical( time_Step , key );
    //Number of attached capture-shrinkage dyneins 
    this->print_dynein_capture_shrinkage_numerical( time_Step , key );
    //Cortical sliding dynein
    this->print_dynein_cortical_sliding_numerical( time_Step , key );
}

//Saves the positions of all microtubule beads and all dyneins
//This should be used very sporadically, since often use would result in very large files
//File name contains the id of the process
void Cell::BIG_PRINT_numerical_micro( unsigned int time_Step , unsigned int key )
{
    //Saves the configuration of the entire cytoskeleton
    this->print_microtubules_numerical( time_Step , key );
    this->print_dynein_numerical( time_Step , key );
}





//This saves the front center of the IS
void Cell::print_center_numerical_parameters( unsigned int key )
{
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "parameters_" + number + ".txt";

      Vector3d center = this->IS_Capture_Shrinkage.get_center_of_IS_front();
      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<center( 0 )<< " " << center( 1 ) <<" "<<center( 2 )<<endl;  //write to file
      fout.close();


}


//Saves the time into the file 
void Cell::print_time( unsigned int time_Step )
{
	FILE *out;
	out = fopen( "./picturesVideos/distance_results/time.txt" , "a");
        double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
        fprintf(out,"%15.12f\n", time );
        fclose( out );
}

//Saves the center of the MTOC into the file
void Cell::print_center_MTOC()
{
	  FILE *out;
	  out = fopen( "./picturesVideos/distance_results/MTOC_center.txt" , "a");
	  Vector3d center = this->MTOC.get_center();
	  fprintf(out,"%15.12f %15.12f %15.12f\n", center( 0 ) , center( 1 ) , center( 2 ) );
	  fclose( out );
}


//Saves the center of the MTOC into the file containing the id key of the process
//File name contains the id of the process
void Cell::print_center_MTOC_numerical( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "MTOC_center_" + number + ".txt";


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<time<<" "<<center( 0 )<< " " << center( 1 ) <<" "<<center( 2 )<<endl;  //write to file
      fout.close();
}


//Saves the number of attached capture-shrinkage dyneins into a file containing the name of the process
//File name contains the id of the process
void Cell::print_dynein_capture_shrinkage_numerical( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "IS_CAP_SH_" + number + ".txt";
      unsigned int number_of_dynein = this->get_number_of_dynein_motors_in_micro_capture();


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<number_of_dynein<<endl;  //write to file
      fout.close();
}




//Saves the number of cortical sliding dyneins acting on MTs
//File name contains the id of the process
void Cell::print_dynein_cortical_sliding_numerical( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "IS_CORTICAL_SLIDING_" + number + ".txt";
      unsigned int number_of_dynein = this->get_number_of_dynein_motors_cortical_sliding_micro();


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<number_of_dynein<<endl; 
      fout.close();
}











//Saves the positions of all dyneins - attached ot unattached
//It saves the dynein to the file whose name contains the i of the process
void Cell::print_dynein_numerical( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      std::string path = "./picturesVideos/numerical_results/";
      std::string number_string = std::to_string( key );
      std::string time_string = std::to_string( time_Step );
      std::string name_of_text_file = path + "dynein_" + time_string + "_"+ number_string + ".txt";

      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      std::map< unsigned int , std::vector<Vector3d> > mapa = this->density_surface_dynein.get_map();
      //Unattached cortical sliding dynein
      for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa.begin(); it != mapa.end(); ++it)
      {

          std::vector<Vector3d> vector_tmp = it->second;
          for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
          {
                Vector3d vect = vector_tmp.at( i );
      		fout<<vect( 0 )<<" "<< vect( 1 ) <<" "<<vect( 2 )<<endl;  //write to file
          }
      }
      fout<<"END"<<endl;

      //Attached cortical sliding dynein
      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( this->array_Of_Microtubules[ micro ].get_dynein_index() == 9 )
          {
                std::vector< std::pair < Vector3d ,double  > > vector_dynein_abscissa = this->array_Of_Microtubules[ micro ].get_Dynein_abscissa();
                if( vector_dynein_abscissa.size() < 1 )
                {
                    continue;
                }
                else
                {
                      for( unsigned int point = 0 ; point < vector_dynein_abscissa.size() ; point ++ )
                      {
                          std::pair < Vector3d ,double  > pair_tmp = vector_dynein_abscissa.at( point );
                          Vector3d position = std::get< 0 >( pair_tmp );
                          fout<<position( 0 )<<" "<< position( 1 ) <<" "<<position( 2 )<<endl;
                      }
                }
          }
      }
      
      fout<<"END"<<endl;
      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( this->array_Of_Microtubules[ micro ].get_dynein_index() == 20 )
          {
                std::vector< std::pair < Vector3d ,double  > > vector_dynein_abscissa = this->array_Of_Microtubules[ micro ].get_Dynein_abscissa();
                if( vector_dynein_abscissa.size() < 1 )
                {
                    continue;
                }
                else
                {
                      for( unsigned int point = 0 ; point < vector_dynein_abscissa.size() ; point ++ )
                      {
                          std::pair < Vector3d ,double  > pair_tmp = vector_dynein_abscissa.at( point );
                          Vector3d position = std::get< 0 >( pair_tmp );
                          fout<<position( 0 )<<" "<< position( 1 ) <<" "<<position( 2 )<<endl;
                      }
                }
          }
      }
      fout.close();
}







//This function saves the entire configuration of all microtubules in cytoskeleton
void Cell::print_microtubules_numerical( unsigned int time_Step , unsigned int key )
{
      //name of the file
      std::string path = "./picturesVideos/numerical_results/";
      std::string number_string = std::to_string( key );
      std::string time_string = std::to_string( time_Step );
      std::string name_of_text_file = path + "micro_" + time_string + "_"+ number_string + ".txt";
      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );

      for( unsigned int microtubule_index = 0 ; microtubule_index < this->number_of_microtubules ; microtubule_index ++ )
      {
          for( unsigned int bead_index = 0 ; bead_index < this->array_Of_Microtubules[ microtubule_index ].getNumberOfPoints() ; bead_index ++ )
          {
                Vector3d micro_point = this->array_Of_Microtubules[ microtubule_index ].getPoint( bead_index );
                fout<<micro_point( 0 )<<" "<< micro_point( 1 ) <<" "<<micro_point( 2 )<<endl;
          }
          fout<<this->array_Of_Microtubules[ microtubule_index ].get_dynein_index() + this->array_Of_Microtubules[ microtubule_index ].get_growth_index( )<<endl;
	  fout<<"END"<<endl;
      }

      fout.close();
}






//The force from the wall acting on every bead of the microtubule
void Cell::force_whole_microtubule_cell_wall( MatrixXd &force , unsigned int index )
{
	if( index >= this->number_of_microtubules )
	{
		cout<<"index >= this->number_of_microtubules in Cell::force_whole_microtubule( MatrixXd &force , unsigned int index )"<<endl;
		cout<<"index = "<<index<<endl;
		throw("");
	}
	Microtubule temporary_microtubule = this->getMicrotubule( index );

	if( force.rows() != 3 * temporary_microtubule.getNumberOfPoints()  )
	{
		cout<<"force.rows() != 3 * temporary_microtubule.getNumberOfPoints() in force_whole_microtubule"<<endl;
		cout<<"index = "<<index<<endl;
		throw("");
	}


	for( unsigned int bead_index = 0 ; bead_index < temporary_microtubule.getNumberOfPoints() ; bead_index ++ )
	{
		Vector3d position = temporary_microtubule.getPoint( bead_index );
        	Vector3d force_one_bead = this->force_position_cell_wall_two_elipsoid( position );
       		for( unsigned int index_dimension = 0 ; index_dimension < 3 ; index_dimension ++ )
		{
			force( 3 * bead_index + index_dimension , 0 ) = force_one_bead( index_dimension );
		}
	}
	if( this->array_Of_Microtubules[ index ].getNumberOfPoints() == 20 )
    	{
        	for( unsigned int index_dimension = 0 ; index_dimension < 3 ; index_dimension ++ )
		{
			force( 3 * ( this->array_Of_Microtubules[ index ].getNumberOfPoints() - 1 ) + index_dimension , 0 ) = 0;
		}
    	}
}












void Cell::timeDevelopment_numerical_results_2( double Time_0 , double Time_1 , unsigned int key )
{

    	double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++

	unsigned int print_integer = 0;

	this->print_center_numerical_parameters( key );
	unsigned int printing_counter = 0;

        //this->BIG_PRINT_numerical_micro( printing_counter , key );
	//printing_counter = printing_counter + 1;
	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{

        	double time = step * sim_of_Cell::time_Step;
        	if( step % 1000 == 0)
		{
            		//cout<<"step % 1000 == 0"<<endl;
            		cout<<"key = "<<key<<endl;
            		cout<<"time = "<<time<<endl;
        	}


      		if( step % 50 == 0 )
		{
			this->BIG_PRINT_numerical( step , key );
        	}

		//ISCorticalSl tmp_1 = this->get_IS_cortical_sliding_first();
		//unsigned int number_caught_dynein = get_number_of_dynein_motors_cortical_sliding_micro();
		//Surface surface( "Cortical_Sliding" , tmp_1 , number_caught_dynein );
		//this->set_density_surface_dynein( surface );



      		if( ( step %  200 ) == 0 ) //14285
		{


/*
			if ( sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage > 0.00001 )
			{
				this->Capture_Shrinkage_dynein.change_part_of_IS_points( this->IS_Capture_Shrinkage );
			}
*/
			if ( sim_of_Cell::density_of_dynein_motor_surface < 0.00001 )
			{	//I fear double imprecisions
				/*
				ISCorticalSl tmp_1 = this->get_IS_cortical_sliding_first();
				unsigned int number_caught_dynein = get_number_of_dynein_motors_cortical_sliding_micro();
				Surface surface( "Cortical_Sliding" , tmp_1 , number_caught_dynein );
				this->set_density_surface_dynein( surface );
				*/
				this->density_surface_dynein.change_part_of_IS_points( this->IS_Cortic_Sl );
			}
			else
			{
				this->density_surface_dynein.change_part_of_points_surface_and_IS( this->IS_Cortic_Sl );
			}
        	}



      		if( ( step %  7142 ) == 0 )
		{
			this->BIG_PRINT_numerical_micro( printing_counter , key );
			printing_counter = printing_counter + 1;
        	}

		//FUNCTIONS
        	//////////////////////////////////////
        	this->MidStep_3(); //

        	//////////////////////////////////////

        	this->stepping_and_detachment_of_all_microtubule_projection_real_dynein();


		//tohle jen kontroluje, jestli je micro chycena v IS
        	this->check_caught_micro_IS_with_real_dynein_all_micro( );
		//tohle pridava dalsi pary do mikrotubule pro IS capt shrinkage
        	this->microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro( );

		this->catch_pair_abscissa_real_dynein();


		//this->control_length_of_micro_IS();
		//this->control_length_of_micro_IS_2();
	}

}














Cell::~Cell()
{

	delete[] array_Of_Microtubules;
	array_Of_Microtubules = NULL;
	// TODO Auto-generated destructor stub
}
