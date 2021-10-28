/*
 * Cell.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: hornak
 */

#include "Cell_two_IS.h"
#include <tuple>



//This is constructor for an unconstrained cytoskeleton in 3D space. Not used for the published results
Cell_two_IS::Cell_two_IS(  )
{
	this->a_axis = Cell_parametres::A_AXIS;
	this->b_axis = Cell_parametres::B_AXIS;
        this->number_of_microtubules_extra = 0;
	this->nucleus = Nucleus();

	//initializing number of microtubules
	this->number_of_microtubules = 16;
	//create array of microtubules
	this->array_Of_Microtubules = new Microtubule[ this->number_of_microtubules ];

	//this is unconstrained space - radius_of_Cell has no meaning
	//IS_position is in x z plane - angle is determined in sim_of_Cell::IS_angle
	//In unconstrained cell IS has no meaning

	this->IS_Capture_Shrinkage = IS_Capt_Shrinkage();
	//this->IS_Cortical_Sliding = ISCorticalSliding();
        this->IS_Cortic_Sl = ISCorticalSl();

        this->MTOC = MTOC2( this->number_of_microtubules );
        //this->abstract_MTOC = ideal_MTOC();
	this->cell_id = 0;
}






//Constructor of the cell with the two IS
//Number of microtubules is given by the file simulationofCell.h - 100
//The variable number_Of_extra_Microtubules would be used to study the asymetrical cytoskeleton.
//The results of the study were not published and number_Of_extra_Microtubules = 100
Cell_two_IS::Cell_two_IS( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant ): Cell( number_Of_Microtubules , number_Of_extra_Microtubules , servant )
{


	//The constructor Cell( number_Of_Microtubules , number_Of_extra_Microtubules , servant ) constructs the cytoskeleton and the first IS
	//This constructor just sets the second IS and simulates dynamic instability
	this->cell_id = servant;
	this->number_of_IS = 2;

	//HERE EVERYTHING ABOUT IS - MICROTUBULE CATCHING WILL BE SET
	//The position of capture shrinkage is set by the orientation
	double azimutal_angle_2 = IS_Capture_shrinkage_param::azimutal_angle_2;
	double planar_component_2 = sin( IS_Capture_shrinkage_param::polar_angle_2 );
	double vertical_component_2 = cos( IS_Capture_shrinkage_param::polar_angle_2 );
	double planar_component_x_2 = planar_component_2 * cos( azimutal_angle_2 );
	double planar_component_y_2 = planar_component_2 * sin( azimutal_angle_2 );

	Vector3d orientation_2( 0.0 , 0.0 , 0.0 );
	orientation_2( 0 ) = planar_component_x_2;
	orientation_2( 1 ) = planar_component_y_2;
	orientation_2( 2 ) = vertical_component_2;
	orientation_2 = orientation_2 / orientation_2.norm();


	//////////////////////////////////////////////////////////CAPTURE-SHRINKAGE 2////////////////////////////////////////////////////////
	
	double radius_of_IS_argument_2 = IS_Capture_shrinkage_param::radius_2;
	Vector3d  center_of_IS_cap_sh_fron_center_2 = IS_Capture_shrinkage_param::z_coordinate_front * orientation_2;

	Vector3d  center_of_IS_cap_sh_rear_center_2 = IS_Capture_shrinkage_param::z_coordinate_back * orientation_2;
	this->IS_Capture_Shrinkage_2 = IS_Capt_Shrinkage( center_of_IS_cap_sh_fron_center_2 , radius_of_IS_argument_2 , center_of_IS_cap_sh_rear_center_2 );


	//////////////////////////////////////////////////////////CORTICAL SLIDING 2////////////////////////////////////////////////////////

    	Vector3d front_center_2 = orientation_2 * IS_Cortical_Sl_parameter::front_radius;
    	Vector3d rear_center_2 = orientation_2 * IS_Cortical_Sl_parameter::rear_radius;

    	double radius_outer_2 = IS_Cortical_Sl_parameter::radius_2;
    	double radius_inner_2 = IS_Cortical_Sl_parameter::radius_inner_2;
    	this->IS_Cortic_Sl_2 = ISCorticalSl( front_center_2 , rear_center_2 , radius_outer_2 , radius_inner_2 );


	//////////////////////////////////////////////////////////DENSITY CAPTURE-SHRINKAGE////////////////////////////////////////////////////////
    	double capture_shrinkage_density = sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage;
    	double capture_shrinkage_density_2 = sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage_2;

    	Surface Surface_capture_shrinkage( capture_shrinkage_density, 1e-6 , this->IS_Capture_Shrinkage );
    	Surface Surface_capture_shrinkage_second( capture_shrinkage_density_2, 1e-6 , this->IS_Capture_Shrinkage_2 );

	//Setting the map holding the points
    	this->set_density_IS_capture_shrinkage( Surface_capture_shrinkage );
    	this->set_density_IS_capture_shrinkage_second( Surface_capture_shrinkage_second );

	//////////////////////////////////////////////////////////DENSITIES CORTICAL SLIDING////////////////////////////////////////////////////////

	//Setting the map holding the points
	 double surface_density =  sim_of_Cell::density_of_dynein_motor_surface;
	 double density_IS_1 = sim_of_Cell::density_of_dynein_motor_Cortical_Sliding;
	 double density_IS_2 = sim_of_Cell::density_of_dynein_motor_Cortical_Sliding_2;
	 Surface first_Surface(  surface_density , sim_of_Cell::surface_width , density_IS_1 , this->get_IS_cortical_sliding_first() , density_IS_2 , this->get_IS_cortical_sliding_second() );
    	 this->set_density_surface_dynein( first_Surface );


	//Dynamic instability ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if( dynamic_instability::instability_switch == false )
	{
		return;
	}


   	unsigned int side = 4;
	unsigned int micro_per_side = this->number_of_microtubules / side;
	unsigned int polygon_per_side;
	   //polygon_per_side gives me the number of polygon at every end of the centrioles

	if( micro_per_side % MTOCparam::micro_in_polygon == 0 )
	{
	    polygon_per_side = micro_per_side / MTOCparam::micro_in_polygon;
	}
	else
	{
	    polygon_per_side = micro_per_side / MTOCparam::micro_in_polygon + 1;
	}
        unsigned int unfinished_polygon_number = micro_per_side % MTOCparam::micro_in_polygon;
	unsigned int counter_micro = 0;
	unsigned int polygon_counter = 0;




	unsigned int Number_of_MT = this->get_microtubule_number();
	MatrixXd lenghts_of_micro_after_DI = MatrixXd::Zero( Number_of_MT , 2 );
	
	//The lengths of microtubules get from the probability distribution
	if( dynamic_instability::decision == 0 )
	{
		lenghts_of_micro_after_DI = this->dynamic_instability_MT_grow(  );
	}
	else if( dynamic_instability::decision == 1 )
	{
		lenghts_of_micro_after_DI = this->dynamic_instability_MT_grow_linear(  );
	}
	else if( dynamic_instability::decision == 2 )
	{
		//This is probability density from the publication
		lenghts_of_micro_after_DI = this->dynamic_instability_MT_grow_August(  );
	}

	if( ( lenghts_of_micro_after_DI.rows(  )  != this->number_of_microtubules ) || ( lenghts_of_micro_after_DI.cols(  )  != 2 ) )
	{
		cout<<"( lenghts_of_micro_after_DI.rows(  )  != this->number_of_microtubules ) || ( lenghts_of_micro_after_DI.cols(  )  != 2 ) "<<endl;
		cout<<"Cell_two_IS::Cell_two_IS( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant ): Cell( number_Of_Microtubules , number_Of_extra_Microtubules , servant )"<<endl;
		cout<<"cols = "<<lenghts_of_micro_after_DI.cols(  )<<endl;
		cout<<"rows = "<<lenghts_of_micro_after_DI.rows(  )<<endl;
		throw("");
	} 

	//Reconstruction of the cytoskeleton with the dynamic instability
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
			   //cout<<"----------------------------------------------------------------------------microtubule = "<<microtubule<<endl;
		           unsigned int point_number = strana * MTOCparam::micro_in_polygon  + microtubule + 1;
			   unsigned int opposite_point_number = this->MTOC.get_random_point_on_opposite_side_without_bias( point_number );
		           Vector3d MTOC_Point = this->MTOC.get_point( point_number );
		           Vector3d opposite_MTOC_Point = this->MTOC.get_point( opposite_point_number );


		           this->abstract_MTOC.set_mtoc_point(  counter_micro , MTOC_Point  );
		           Vector3d orientation = MTOC_Point - opposite_MTOC_Point;
		           orientation = orientation / orientation.norm();


		           Vector3d first_Point = opposite_MTOC_Point;
		           Vector3d second_Point = MTOC_Point;

			   double Lenght_of_micro = lenghts_of_micro_after_DI( counter_micro , 0 );
			   if( Lenght_of_micro < dynamic_instability::coefficient_4 )
			   {
				cout<<" Lenght_of_micro < dynamic_instability::coefficient_4"<<endl;
				cout<<"Lenght_of_micro = "<<Lenght_of_micro<<endl;
				throw("");
			   }		

			   //Index determining, whether the microtubule grows, or shrinks
			   double growing_index = lenghts_of_micro_after_DI( counter_micro , 1 );
			
			   unsigned int number_of_points_tmp = ( unsigned int ) ( floor( Lenght_of_micro / sim_of_Cell::resting_distance ) );	
			   unsigned int number_of_points_tmp_2 = number_of_points_tmp;

			   if( number_of_points_tmp == 0 )
			   {
				number_of_points_tmp = 1;
			   } 	 	 		
			   unsigned int number_of_points = number_of_points_tmp + 2;
			   //The microtubule is created
		           Microtubule tmp( first_Point , second_Point , orientation , microtubule , polygon_counter , strana , point_number , opposite_point_number , number_of_points );//true
			   tmp.set_Cell_ID( this->cell_id );	
	
			   if( number_of_points_tmp_2 >= 1 )
			   {
				double additional_lenght = Lenght_of_micro - sim_of_Cell::resting_distance * ( double ) number_of_points_tmp_2;
				tmp.prolong_Last_tangent( additional_lenght );
				//Resizing of the microtubule
				tmp.control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0();
			   }
			   else if( number_of_points_tmp_2 == 0 )
			   {				
				double shortening = sim_of_Cell::resting_distance - Lenght_of_micro;
				tmp.shorten_Last_tangent( shortening);
			   } 	

			   unsigned int growing_index_tmp = 0;
			   if( abs( growing_index - 1.0 ) < 1.0e-6 )
			   {
				growing_index_tmp = 1;
			   }		
			   else if( abs( growing_index - 2.0 ) < 1.0e-6 )
			   {
				growing_index_tmp = 2;
			   }		
			   else
			   {
			    cout<<"Cell_two_IS( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant ): Cell( number_Of_Microtubules , number_Of_extra_Microtubules , servant )"<<endl;
			    cout<<"growing_index =  "<<growing_index<<endl;
			    cout<<"ERROR_ID = 6749684564561"<<endl;
			    throw("");
			   }		 	

			   tmp.set_growth_index( growing_index_tmp );
			   tmp.set_lenght_of_micro_0_outside_MTOC( tmp.get_lenght_of_microtubule_outside_MTOC() );

		           this->array_Of_Microtubules[ counter_micro ] = tmp;
                      	   if( this->array_Of_Microtubules[ counter_micro ].get_last_Tangent( ).norm() < dynamic_instability::coefficient_4 )
			   {

			   }

			   //Abstract MTOC is not used - the object for later simulations
		           this->abstract_MTOC.set_orientation( counter_micro , this->array_Of_Microtubules[ counter_micro ].getPoint( 1 ) );
		           Vector3d mtoc_point = this->MTOC.get_center() - orientation * sim_of_Cell::resting_distance;
		           this->abstract_MTOC.set_point( counter_micro , mtoc_point );
		           counter_micro = counter_micro + 1;
		    }
		    polygon_counter = polygon_counter + 1;
		}

	 }










}





//Basic setter
//This just takes a created Surface object and give it to the variable
void Cell_two_IS::set_density_IS_capture_shrinkage_second( Surface& density )
{
	if( this->number_of_IS != 2 )
	{
		cout<<" this->number_of_IS != 2 "<<endl;
		cout<<" Cell::set_density_IS_capture_shrinkage_second() "<<endl;
		cout<<"ERROR_ID = 9494646164"<<endl;
		throw("");
	}
	this->Capture_Shrinkage_dynein_2 = density;
}


//Gets the second capture-shrinkage IS
IS_Capt_Shrinkage Cell_two_IS::get_IS_capture_shrinkage_second()
{
	if( this->number_of_IS != 2 )
	{
		cout<<" this->number_of_IS != 2 "<<endl;
		cout<<" Cell::get_IS_capture_shrinkage_second() "<<endl;
		cout<<"ERROR_ID = 49861561"<<endl;
		throw("");
	}
	return this->IS_Capture_Shrinkage_2;

}







//This function makes files for the generation of the video
//Basically, this functions takes every point and all the geometrical objects and handles them for the Python for the preparation of Povray files
//They are saved in the way that Povray can easily display the whole cell
void Cell_two_IS::print_Cell_two_IS( unsigned int index )
{
          char name_of_text_file [55];
	  sprintf ( name_of_text_file , "picturesVideos/textFiles/microtubule_%d.txt", index );
	  FILE *out;
	  out = fopen( name_of_text_file , "w");


          //The printing of the microtubules
          //Unattached, attached MTs are distinguished by the index
          //This will later allow Povray to use different colors
	  for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
	  {
		if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 0 )
         	{
	                Vector3d center = this->MTOC.get_center();
	                double parametr = 0.0;
                        fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  center( 0 )  * 1e6, center( 1 )  * 1e6, center( 2 )   * 1e6);			
			
			//Beads of the microtubules are written to the file
			//Increasing parameter before 3D point is required by Povray 
			//The coordinates are written in metres, since Povray dislikes 1e-6
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
                	//I generate one extra point, since Povray splines never prints the last point
                	fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                	//The following integer is used to distinguish between unattached and attached microtubules so that the Povray can use diferent colors
                	//The case of unattached microtubule is special, since the growing and shrinking MTs are printed with different colors
                	fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() + this->array_Of_Microtubules[ microtubule ].get_growth_index( ) );
                     	fprintf(out,"END\n" );
	  }
	 //The same is done for different Microtubules 
         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
         {
		    //continue;
                    Vector3d center = this->MTOC.get_center();
                    double parametr = 0.0;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  center( 0 )  * 1e6, center( 1 )  * 1e6, center( 2 )   * 1e6);

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
                    fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() + this->array_Of_Microtubules[ microtubule ].get_growth_index( ) );
                    fprintf(out,"END\n" );

         }
         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
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
         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 )
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
         else
	 {
	 	//Just checks whether a microtubule has a different, unexpected index
	 	//If this is the case, the simulation is stoped by error and the description
		cout<<endl;
		cout<<"this->array_Of_Microtubules[ microtubule ].get_dynein_index() = "<<this->array_Of_Microtubules[ microtubule ].get_dynein_index() <<endl;
		cout<<endl;
		cout<<endl;
		throw("");	
	 }

	  }
	  fclose( out );



          //MTOC printing
	  char name_of_MTOC_file [55];
	  sprintf ( name_of_MTOC_file , "picturesVideos/textFiles/MTOC_%d.txt", index );
	  FILE *out_MTOC;
	  out_MTOC = fopen( name_of_MTOC_file , "w");
	  //First, the center is printed
	  //Again, the coordinates are multiplied by 1e6 due to Povray 
          Vector3d MTOC_point = this->MTOC.get_center();
          fprintf( out_MTOC , "%10.5f,%10.5lf,%10.5lf\n",  MTOC_point( 0 ) * 1e6 , MTOC_point( 1 ) * 1e6 , MTOC_point( 2 ) * 1e6 );
          fprintf( out_MTOC , "END\n" );
	  //All the points are printed
	  for( unsigned int poly_number = 1 ; poly_number <= this->MTOC.get_number_of_points() ; poly_number ++ )
	  {
                  Vector3d MTOC_point =   this->MTOC.get_point( poly_number );
                  //Since the MTOC is not printed using the slides, parameter is not required
		  fprintf( out_MTOC , "%10.5f,%10.5lf,%10.5lf\n",  MTOC_point( 0 ) * 1e6 , MTOC_point( 1 ) * 1e6 , MTOC_point( 2 ) * 1e6 );
		  //The END is delimiter for the Python
		  fprintf( out_MTOC , "END\n" );
	  }
	  fclose( out_MTOC );


	//This file contains only the center of the MTOC 
      	char MTOC_center [55];
      	sprintf ( MTOC_center , "picturesVideos/textFiles/MTOC_center_%d.txt", index );
      	FILE *out_MTOC_center;
      	out_MTOC_center = fopen( MTOC_center , "w");
      	//I get the center, print it, close the file
      	Vector3d center_of_MTOC = this->MTOC.get_center();
      	fprintf( out_MTOC_center , "%10.5f,%10.5lf,%10.5lf\n",  center_of_MTOC( 0 ) * 1e6 , center_of_MTOC( 1 ) * 1e6 , center_of_MTOC( 2 ) * 1e6 );
      	fclose( out_MTOC_center );


	  //This prints the coordinates of the points, where the plus-end of the microtubule depolimerizes
	  char name_of_IS_cathing_file [55];
	  sprintf ( name_of_IS_cathing_file , "picturesVideos/textFiles/IS_catching_%d.txt", index );
	  FILE *out_IS_catching;
	  out_IS_catching = fopen( name_of_IS_cathing_file , "w");
	  for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
	  {
	        //It checks whether the microtubule is attached to the capture-shrinkage dyneins in the first, or the second IS
	        //First, I have to check whether the microtubules are attached to capture-shrinkage dyneins 
	   	if( ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 ) || ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 ) )
		{
            		Vector3d IS_catching = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();
			fprintf( out_IS_catching , "%10.5f,%10.5lf,%10.5lf\n",  IS_catching( 0 ) * 1e6 , IS_catching( 1 ) * 1e6 , IS_catching( 2 ) * 1e6 );
                	fprintf( out_IS_catching , "END\n" );
            	}

          }
	  fclose( out_IS_catching );



	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //The geometry of the IS are printed
	  //This prints IS capture-shrinkage geometry
	  //The center of the front and the rear side of the IS are printed with the radius of the IS
	  //The Povray prints cylinder from the two points and the radius
	  char IS_capture_shrinkage [55];
	  sprintf ( IS_capture_shrinkage , "picturesVideos/textFiles/IS_capture_shrinkage_%d.txt", index );
	  FILE *out_IS_capture_shrinkage;
	  out_IS_capture_shrinkage = fopen( IS_capture_shrinkage , "w");
	  
	 //I read the centers and the radius of the IS 1 capture-shrinkage
	 //Center of the front side
	  Vector3d IS_capture_shrinkage_center = this->IS_Capture_Shrinkage.get_center_of_IS_front();
	  double radius_IS_Capture_Shrinkage = this->IS_Capture_Shrinkage.get_radius_of_IS();
	  //Center of the rear side of the IS
	  Vector3d IS_capture_shrinkage_rear_point = this->IS_Capture_Shrinkage.get_center_of_IS_rear();


	  //They are printed in the way that the Python can easily read them
	  fprintf( out_IS_capture_shrinkage , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_center( 0 ) * 1e6 , IS_capture_shrinkage_center( 1 ) * 1e6 , IS_capture_shrinkage_center( 2 ) * 1e6 );
	  fprintf( out_IS_capture_shrinkage , "END\n" );
	  fprintf( out_IS_capture_shrinkage , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_rear_point( 0 ) * 1e6 , IS_capture_shrinkage_rear_point( 1 ) * 1e6 , IS_capture_shrinkage_rear_point( 2 )  * 1e6 );
	  fprintf( out_IS_capture_shrinkage , "END\n" );
	  fprintf( out_IS_capture_shrinkage , "%10.5f\n",  radius_IS_Capture_Shrinkage * 1e6 );
          fclose( out_IS_capture_shrinkage );

	  //Printing of the cortical sliding IS in the same way as the capture-shrinkage IS	
          char IS_cortical_sl[55];
	  sprintf ( IS_cortical_sl , "picturesVideos/textFiles/IS_cortical_sl_%d.txt", index );
	  FILE *out_IS_cortical_sl;
	  out_IS_cortical_sl = fopen( IS_cortical_sl , "w");

	  Vector3d front_center1 = this->IS_Cortic_Sl.get_center_front_of_IS();
	  Vector3d rear_center1 = this->IS_Cortic_Sl.get_center_rear_of_IS();

	  Vector3d point_on_plane = this->IS_Cortic_Sl.get_point_on_plane();
	  Vector3d axis = this->IS_Cortic_Sl.get_axis();
          Vector3d tmp_rear_point = point_on_plane + ( - 1.0 ) * axis * 2e-6;
	  double radius = this->IS_Cortic_Sl.get_radius();
          fprintf( out_IS_cortical_sl , "%10.5f,%10.5lf,%10.5lf\n",  front_center1( 0 ) * 1e6 , front_center1( 1 ) * 1e6 , front_center1( 2 ) * 1e6 );
	  fprintf( out_IS_cortical_sl , "END\n" );
	  fprintf( out_IS_cortical_sl , "%10.5f,%10.5lf,%10.5lf\n",  rear_center1( 0 ) * 1e6 , rear_center1( 1 ) * 1e6 , rear_center1( 2 )  * 1e6 );
	  fprintf( out_IS_cortical_sl , "END\n" );
	  fprintf( out_IS_cortical_sl , "%10.5f\n",  radius * 1e6 );
	  fclose( out_IS_cortical_sl );





      /////////////////////////////////////////////////////////////////////////////////////////////
      //This prints the dynein distributed in both cortical sliding IS
      //The dynein is saved in the std maps
      //To print the dyneins, we get the map and iterate through it
           char Dynein_surface_randomly_distributed [155];
           sprintf ( Dynein_surface_randomly_distributed , "picturesVideos/textFiles/Dynein_surface_randomly_distributed_%d.txt", index );
           FILE *out_Dynein_surface_randomly_distributed;
           out_Dynein_surface_randomly_distributed = fopen( Dynein_surface_randomly_distributed , "w");
           std::map< unsigned int , std::vector<Vector3d> > mapa = this->density_surface_dynein.get_map();
           for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa.begin(); it != mapa.end(); ++it)
           {
		//We get the vector of 3D points belonging to the place in the map
                std::vector<Vector3d> vector_tmp = it->second;
                //For loop goes through the vector and prints the points
                for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
                {
                    Vector3d vect = vector_tmp.at( i );
                    fprintf( out_Dynein_surface_randomly_distributed , "%10.5f,%10.5lf,%10.5lf\n",  vect( 0 ) * 1e6 , vect( 1 ) * 1e6 , vect( 2 ) * 1e6 );
                    fprintf( out_Dynein_surface_randomly_distributed , "END\n" );
                }
           }
           fclose( out_Dynein_surface_randomly_distributed );
           
           
          //This prints the dynein distributed in IS capture shrinkage
          //This takes the std:map holding the coordinates of the unattached dyneins represented only by 3D points  
           char Dynein_IS_capture[155];
           sprintf ( Dynein_IS_capture , "picturesVideos/textFiles/Dynein_IS_capture_%d.txt", index );
           FILE *out_Dynein_IS_capture;
           out_Dynein_IS_capture = fopen( Dynein_IS_capture , "w");
           std::map< unsigned int , std::vector<Vector3d> > mapa_capt = this->Capture_Shrinkage_dynein.get_map();
           for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa_capt.begin(); it != mapa_capt.end(); ++it)
           {
                std::vector<Vector3d> vector_tmp = it->second;
                //it->second holds the std vector of the 3D points
                //The cyckle goes through the vector
                for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
                {
                    Vector3d vect = vector_tmp.at( i );
                    fprintf( out_Dynein_IS_capture , "%10.5f,%10.5lf,%10.5lf\n",  vect( 0 ) * 1e6 , vect( 1 ) * 1e6 , vect( 2 ) * 1e6 );
                    fprintf( out_Dynein_IS_capture , "END\n" );
                }
           }
           fclose( out_Dynein_IS_capture );



      // Prints the anchor points of dyneins attached to microtubules
      // Dyneins attached to microtubules are represented by the anchor points and abscissa, which is the lenght between the attachment point and the beginning of the microtubule 
      // For each microtubule, they are stored in the vector of pairs std::pair < Vector3d ,double  >
      char dynein_abscissa [155];
      sprintf ( dynein_abscissa , "picturesVideos/textFiles/dynein_abscissa_%d.txt", index );
      FILE *out_dynein_abscissa;
      out_dynein_abscissa = fopen( dynein_abscissa , "w");
      //The cycle goes through the microtubules to find the ones with dynein
      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 9 ) || ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 20 )  || ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 40 ))
          {
                //We get the vectors of pairs
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
                          //We get the anchor point
                          Vector3d position = std::get< 0 >( pair_tmp );
                          fprintf(out_dynein_abscissa,"%10.5f,%10.5lf,%10.5lf\n",position( 0 ) * 1e6 , position( 1 ) * 1e6 , position( 2 ) * 1e6 );
                          fprintf( out_dynein_abscissa , "END\n" );
                      }
                }
          }
      }
      fclose( out_dynein_abscissa );

      //Print the attachment points cortical sliding
      //Again, we get the std vector of std pairs
      //The position of the attachment point is calculated from the length of the abscissa
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
                    //I take abscissa
                    double abscissa =   std::get<1>( one_pair );
                    unsigned int lower_bead_index =  this->array_Of_Microtubules[ micro ].get_index_according_to_abscissa( abscissa );

                    Vector3d lower_point = this->array_Of_Microtubules[ micro ].getPoint( lower_bead_index );
                    Vector3d tangent = this->array_Of_Microtubules[ micro ].getTangent2( lower_bead_index );


                    double abscissa_minus_lower_bead = abscissa - ( ( double ) lower_bead_index ) * this->array_Of_Microtubules[ micro ].getRestDist();
                    //The position of the attachment point is calculated
                    Vector3d p_of_att = lower_point + ( abscissa_minus_lower_bead / tangent.norm()  ) * tangent; // this->getRestDist()
                   
                    
                    fprintf(out_dynein_abscissa_attachment,"%10.5f,%10.5lf,%10.5lf\n",p_of_att( 0 )*1e6,p_of_att( 1 )*1e6,p_of_att( 2 )*1e6 );
                    fprintf( out_dynein_abscissa_attachment , "END\n" );

                }
          }

      }
      fclose( out_dynein_abscissa_attachment );
      
      //Here the capture-shrinkage nd cortical sliding IS2 are printed in the same way as IS1, see the comments above
      if( this->number_of_IS == 2 )
      {
	 	//this concerns IS geometry and dyneins

	 	//IS capture shrinkage GEOMETRY
	 	
	  	char IS_capture_shrinkage_2 [55];
	  	sprintf ( IS_capture_shrinkage_2 , "picturesVideos/textFiles/IS_capture_shrinkage_2_%d.txt", index );
	  	FILE *out_IS_capture_shrinkage_2;
	  	out_IS_capture_shrinkage_2 = fopen( IS_capture_shrinkage_2 , "w");
	  	Vector3d IS_capture_shrinkage_center_2 = this->IS_Capture_Shrinkage_2.get_center_of_IS_front();
	  	double radius_IS_Capture_Shrinkage_2 = this->IS_Capture_Shrinkage_2.get_radius_of_IS();
	  	Vector3d IS_capture_shrinkage_rear_point_2 = this->IS_Capture_Shrinkage_2.get_center_of_IS_rear();


		fprintf( out_IS_capture_shrinkage_2 , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_center_2( 0 ) * 1e6 , IS_capture_shrinkage_center_2( 1 ) * 1e6 , IS_capture_shrinkage_center_2( 2 ) * 1e6 );
		fprintf( out_IS_capture_shrinkage_2 , "END\n" );
		fprintf( out_IS_capture_shrinkage_2 , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_rear_point_2( 0 ) * 1e6 , IS_capture_shrinkage_rear_point_2( 1 ) * 1e6 , IS_capture_shrinkage_rear_point_2( 2 )  * 1e6 );
		fprintf( out_IS_capture_shrinkage_2 , "END\n" );
		fprintf( out_IS_capture_shrinkage_2 , "%10.5f\n",  radius_IS_Capture_Shrinkage_2 * 1e6 );
		fclose( out_IS_capture_shrinkage_2 );
		
		
	  //!!!!!!!!!!!!!!!!!!!! CORTICAL SLIDING GEOMETRY

          	char IS_cortical_sl_2[55];
	  	sprintf ( IS_cortical_sl_2 , "picturesVideos/textFiles/IS_cortical_sl_2_%d.txt", index );
	  	FILE *out_IS_cortical_sl_2;
	  	out_IS_cortical_sl_2 = fopen( IS_cortical_sl_2 , "w");

	  	Vector3d front_center_2 = this->IS_Cortic_Sl_2.get_center_front_of_IS();
	  	Vector3d rear_center_2 = this->IS_Cortic_Sl_2.get_center_rear_of_IS();


		  Vector3d point_on_plane_2 = this->IS_Cortic_Sl_2.get_point_on_plane();
		  Vector3d axis_2 = this->IS_Cortic_Sl_2.get_axis();
	          Vector3d tmp_rear_point_2 = point_on_plane_2 + ( - 1.0 ) * axis_2 * 2e-6;
		  double radius_2 = this->IS_Cortic_Sl_2.get_radius();
	          fprintf( out_IS_cortical_sl_2 , "%10.5f,%10.5lf,%10.5lf\n",  front_center_2( 0 ) * 1e6 , front_center_2( 1 ) * 1e6 , front_center_2( 2 ) * 1e6 );
		  fprintf( out_IS_cortical_sl_2 , "END\n" );
		  fprintf( out_IS_cortical_sl_2 , "%10.5f,%10.5lf,%10.5lf\n",  rear_center_2( 0 ) * 1e6 , rear_center_2( 1 ) * 1e6 , rear_center_2( 2 )  * 1e6 );
		  fprintf( out_IS_cortical_sl_2 , "END\n" );
		  fprintf( out_IS_cortical_sl_2 , "%10.5f\n",  radius_2 * 1e6 );
		  fclose( out_IS_cortical_sl_2 );


          	//This prints the dynein distributed in IS2 capture-shrinkage, see the comments above
           	char Dynein_IS_capture_2[155];
           	sprintf ( Dynein_IS_capture_2 , "picturesVideos/textFiles/Dynein_IS_capture_2_%d.txt", index );
           	FILE *out_Dynein_IS_capture_2;
           	out_Dynein_IS_capture_2 = fopen( Dynein_IS_capture_2 , "w");
           	std::map< unsigned int , std::vector<Vector3d> > mapa_capt_2 = this->Capture_Shrinkage_dynein_2.get_map();
           	for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa_capt_2.begin(); it != mapa_capt_2.end(); ++it)
           	{

                	std::vector<Vector3d> vector_tmp = it->second;
                	for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
                	{
                    	Vector3d vect = vector_tmp.at( i );
                    	fprintf( out_Dynein_IS_capture_2 , "%10.5f,%10.5lf,%10.5lf\n",  vect( 0 ) * 1e6 , vect( 1 ) * 1e6 , vect( 2 ) * 1e6 );
                    	fprintf( out_Dynein_IS_capture_2 , "END\n" );
                	}
           	}

           	fclose( out_Dynein_IS_capture_2 );
      }



	   //The nucleus is printed here
           char nucleus_print[155];
           sprintf ( nucleus_print , "picturesVideos/textFiles/nucleus_print_%d.txt", index );
           FILE *out_nucleus_print;
           out_nucleus_print = fopen( nucleus_print , "w");
           //We print the center of the nucleus
           Vector3d nucleus_center = this->nucleus.get_center();
           fprintf( out_nucleus_print , "%10.5f,%10.5lf,%10.5lf\n",  nucleus_center( 0 ) * 1e6 , nucleus_center( 1 ) * 1e6 , nucleus_center( 2 ) * 1e6 );
           //And the two Axis
           //In simulations, axis are the same, since the nucleus is a sphere
	   fprintf( out_nucleus_print , "%10.5f,%10.5lf\n",  this->nucleus.get_A_Axis() * 1e6 , this->nucleus.get_B_Axis() * 1e6 );
           fprintf( out_nucleus_print , "END\n" );
           fclose( out_nucleus_print );

}





//This function does basically the same as the previous one - prints the entire configuration of the cell and the cytoskeleton
//It is called not for the purposes of the video, but for the purposes of the control
//If something goes wrong in the simulation and NaN is detected, the error terminates the simulation and this function is called to save the configuration
void Cell_two_IS::print_Cell_two_IS_numerical( unsigned int index )
{
          char name_of_text_file [55];
	  sprintf ( name_of_text_file , "picturesVideos/numerical_results/microtubule_%d.txt", index );

	  FILE *out;
	  out = fopen( name_of_text_file , "w");

	  for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
	  {
		if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 0 )
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
                	Vector3d ill_point = point_last + tangent_last ; // / 10.0
                	parametr = parametr + 0.05;
                	fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                	fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() + this->array_Of_Microtubules[ microtubule ].get_growth_index( ) );
                     	fprintf(out,"END\n" );
	  }
         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
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
                    Vector3d tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                    Vector3d point_last = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
                    Vector3d ill_point = point_last + tangent_last / 10.0;
                    parametr = parametr + 0.05;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                    fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() + this->array_Of_Microtubules[ microtubule ].get_growth_index( ) );
                    fprintf(out,"END\n" );
         }

         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
         {
		    //continue;
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
         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 )
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
         else
	 {
		cout<<endl;
		cout<<endl;
		cout<<"this->array_Of_Microtubules[ microtubule ].get_dynein_index() = "<<this->array_Of_Microtubules[ microtubule ].get_dynein_index() <<endl;
		cout<<endl;
		cout<<endl;
		throw("");	
	 }

	  }
	  fclose( out );



      //MTOC printing---------------------------------------------------------------------------------------
	  char name_of_MTOC_file [55];
	  sprintf ( name_of_MTOC_file , "picturesVideos/numerical_results/MTOC_%d.txt", index );
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

      	char MTOC_center [55];
      	sprintf ( MTOC_center , "picturesVideos/numerical_results/MTOC_center_%d.txt", index );
      	FILE *out_MTOC_center;
      	out_MTOC_center = fopen( MTOC_center , "w");
      	Vector3d center_of_MTOC = this->MTOC.get_center();
      	fprintf( out_MTOC_center , "%10.5f,%10.5lf,%10.5lf\n",  center_of_MTOC( 0 ) * 1e6 , center_of_MTOC( 1 ) * 1e6 , center_of_MTOC( 2 ) * 1e6 );
      	fclose( out_MTOC_center );


	  char name_of_IS_cathing_file [55];
	  sprintf ( name_of_IS_cathing_file , "picturesVideos/numerical_results/IS_catching_%d.txt", index );
	  FILE *out_IS_catching;
	  out_IS_catching = fopen( name_of_IS_cathing_file , "w");

	  for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
	  {
	   	if( ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 ) || ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 ) )
		{
            		Vector3d IS_catching = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();
			fprintf( out_IS_catching , "%10.5f,%10.5lf,%10.5lf\n",  IS_catching( 0 ) * 1e6 , IS_catching( 1 ) * 1e6 , IS_catching( 2 ) * 1e6 );
                	fprintf( out_IS_catching , "END\n" );
            	}

          }
	  fclose( out_IS_catching );


	  char IS_capture_shrinkage [55];
	  sprintf ( IS_capture_shrinkage , "picturesVideos/numerical_results/IS_capture_shrinkage_%d.txt", index );
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
	  sprintf ( IS_cortical_sl , "picturesVideos/numerical_results/IS_cortical_sl_%d.txt", index );
	  FILE *out_IS_cortical_sl;
	  out_IS_cortical_sl = fopen( IS_cortical_sl , "w");

	  Vector3d front_center1 = this->IS_Cortic_Sl.get_center_front_of_IS();
	  //cout<<front_center<<endl;
	  Vector3d rear_center1 = this->IS_Cortic_Sl.get_center_rear_of_IS();


	  Vector3d point_on_plane = this->IS_Cortic_Sl.get_point_on_plane();
	  Vector3d axis = this->IS_Cortic_Sl.get_axis();
          Vector3d tmp_rear_point = point_on_plane + ( - 1.0 ) * axis * 2e-6;
	  double radius = this->IS_Cortic_Sl.get_radius();
          fprintf( out_IS_cortical_sl , "%10.5f,%10.5lf,%10.5lf\n",  front_center1( 0 ) * 1e6 , front_center1( 1 ) * 1e6 , front_center1( 2 ) * 1e6 );
	  fprintf( out_IS_cortical_sl , "END\n" );
	  fprintf( out_IS_cortical_sl , "%10.5f,%10.5lf,%10.5lf\n",  rear_center1( 0 ) * 1e6 , rear_center1( 1 ) * 1e6 , rear_center1( 2 )  * 1e6 );
	  fprintf( out_IS_cortical_sl , "END\n" );
	  fprintf( out_IS_cortical_sl , "%10.5f\n",  radius * 1e6 );
	  fclose( out_IS_cortical_sl );



           char Dynein_surface_randomly_distributed [155];
           sprintf ( Dynein_surface_randomly_distributed , "picturesVideos/numerical_results/Dynein_surface_randomly_distributed_%d.txt", index );
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
           sprintf ( Dynein_IS_capture , "picturesVideos/numerical_results/Dynein_IS_capture_%d.txt", index );
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


      //print the attachment points in microtubules
      char dynein_abscissa [155];
      sprintf ( dynein_abscissa , "picturesVideos/numerical_results/dynein_abscissa_%d.txt", index );
      FILE *out_dynein_abscissa;
      out_dynein_abscissa = fopen( dynein_abscissa , "w");
      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 9 ) || ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 20 )  || ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 40 ))
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
      sprintf ( dynein_abscissa_attachment , "picturesVideos/numerical_results/dynein_abscissa_attachment_%d.txt", index );
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

      if( this->number_of_IS == 2 )
      {
	 	//this concerns IS geometry and dyneins

	  	char IS_capture_shrinkage_2 [55];
	  	sprintf ( IS_capture_shrinkage_2 , "picturesVideos/numerical_results/IS_capture_shrinkage_2_%d.txt", index );
	  	FILE *out_IS_capture_shrinkage_2;
	  	out_IS_capture_shrinkage_2 = fopen( IS_capture_shrinkage_2 , "w");
	  	Vector3d IS_capture_shrinkage_center_2 = this->IS_Capture_Shrinkage_2.get_center_of_IS_front();
	  	double radius_IS_Capture_Shrinkage_2 = this->IS_Capture_Shrinkage_2.get_radius_of_IS();
	  	Vector3d IS_capture_shrinkage_rear_point_2 = this->IS_Capture_Shrinkage_2.get_center_of_IS_rear();

		fprintf( out_IS_capture_shrinkage_2 , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_center_2( 0 ) * 1e6 , IS_capture_shrinkage_center_2( 1 ) * 1e6 , IS_capture_shrinkage_center_2( 2 ) * 1e6 );
		fprintf( out_IS_capture_shrinkage_2 , "END\n" );
		fprintf( out_IS_capture_shrinkage_2 , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_rear_point_2( 0 ) * 1e6 , IS_capture_shrinkage_rear_point_2( 1 ) * 1e6 , IS_capture_shrinkage_rear_point_2( 2 )  * 1e6 );
		fprintf( out_IS_capture_shrinkage_2 , "END\n" );
		fprintf( out_IS_capture_shrinkage_2 , "%10.5f\n",  radius_IS_Capture_Shrinkage_2 * 1e6 );
		fclose( out_IS_capture_shrinkage_2 );


          	char IS_cortical_sl_2[55];
	  	sprintf ( IS_cortical_sl_2 , "picturesVideos/numerical_results/IS_cortical_sl_2_%d.txt", index );
	  	FILE *out_IS_cortical_sl_2;
	  	out_IS_cortical_sl_2 = fopen( IS_cortical_sl_2 , "w");

	  	Vector3d front_center_2 = this->IS_Cortic_Sl_2.get_center_front_of_IS();
	  	Vector3d rear_center_2 = this->IS_Cortic_Sl_2.get_center_rear_of_IS();


		  Vector3d point_on_plane_2 = this->IS_Cortic_Sl_2.get_point_on_plane();
		  Vector3d axis_2 = this->IS_Cortic_Sl_2.get_axis();
	          Vector3d tmp_rear_point_2 = point_on_plane_2 + ( - 1.0 ) * axis_2 * 2e-6;
		  double radius_2 = this->IS_Cortic_Sl_2.get_radius();
	          fprintf( out_IS_cortical_sl_2 , "%10.5f,%10.5lf,%10.5lf\n",  front_center_2( 0 ) * 1e6 , front_center_2( 1 ) * 1e6 , front_center_2( 2 ) * 1e6 );
		  fprintf( out_IS_cortical_sl_2 , "END\n" );
		  fprintf( out_IS_cortical_sl_2 , "%10.5f,%10.5lf,%10.5lf\n",  rear_center_2( 0 ) * 1e6 , rear_center_2( 1 ) * 1e6 , rear_center_2( 2 )  * 1e6 );
		  fprintf( out_IS_cortical_sl_2 , "END\n" );
		  fprintf( out_IS_cortical_sl_2 , "%10.5f\n",  radius_2 * 1e6 );
		  fclose( out_IS_cortical_sl_2 );

          	//this prints the dynein distributed in IS capture shrinkage
           	char Dynein_IS_capture_2[155];
           	sprintf ( Dynein_IS_capture_2 , "picturesVideos/numerical_results/Dynein_IS_capture_2_%d.txt", index );
           	FILE *out_Dynein_IS_capture_2;
           	out_Dynein_IS_capture_2 = fopen( Dynein_IS_capture_2 , "w");
           	std::map< unsigned int , std::vector<Vector3d> > mapa_capt_2 = this->Capture_Shrinkage_dynein_2.get_map();
           	for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa_capt_2.begin(); it != mapa_capt_2.end(); ++it)
           	{

                	std::vector<Vector3d> vector_tmp = it->second;
                	for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
                	{
                    	Vector3d vect = vector_tmp.at( i );
                    	fprintf( out_Dynein_IS_capture_2 , "%10.5f,%10.5lf,%10.5lf\n",  vect( 0 ) * 1e6 , vect( 1 ) * 1e6 , vect( 2 ) * 1e6 );
                    	fprintf( out_Dynein_IS_capture_2 , "END\n" );
                	}
           	}

           	fclose( out_Dynein_IS_capture_2 );


      }


           char nucleus_print[155];
           sprintf ( nucleus_print , "picturesVideos/numerical_results/nucleus_print_%d.txt", index );
           FILE *out_nucleus_print;
           out_nucleus_print = fopen( nucleus_print , "w");
           Vector3d nucleus_center = this->nucleus.get_center();
           fprintf( out_nucleus_print , "%10.5f,%10.5lf,%10.5lf\n",  nucleus_center( 0 ) * 1e6 , nucleus_center( 1 ) * 1e6 , nucleus_center( 2 ) * 1e6 );
	   fprintf( out_nucleus_print , "%10.5f,%10.5lf\n",  this->nucleus.get_A_Axis() * 1e6 , this->nucleus.get_B_Axis() * 1e6 );
           fprintf( out_nucleus_print , "END\n" );
           fclose( out_nucleus_print );

}

















//If not a number detected, the simulation is immediately stopped 	
void Cell_two_IS::control_NAN( string arg )
{
    for( unsigned int micro_index = 0 ; micro_index < this->number_of_microtubules ; micro_index ++ )
    {
	for( unsigned int i_tmp = 0 ; i_tmp <  this->array_Of_Microtubules[ micro_index ].getNumberOfPoints() ; i_tmp++ )
	{
		Vector3d point = this->array_Of_Microtubules[ micro_index ].getPoint( i_tmp );
		if(  std::isnan( point.norm() ) == true  )
		{
			cout<<"std::isnan( point.norm() ) == true "<<endl;
			cout<<"micro_index = "<<micro_index<<endl;
			cout<<"i_tmp = "<<i_tmp<<endl;
			cout<<"this->array_Of_Microtubules[ micro_index ].get_dynein_index() = "<<this->array_Of_Microtubules[ micro_index ].get_dynein_index()<<endl;
			cout<<"this->array_Of_Microtubules[ micro_index ].get_growth_index( ) = "<<this->array_Of_Microtubules[ micro_index ].get_growth_index( )<<endl;
			cout<<"arg = "<<arg<<endl;
			cout<<"this->cell_id = "<<this->cell_id<<endl;
			this->print_Cell_two_IS_numerical( 1000 );
			throw("");
			
		}  

	}


    }	
}






//This function resizes the cytoskeleton after every step
//The projection operator project forces on hyperplane so that they are perpendicular to the constraints
//Still, the lenghts of the microtubules can change due to numerical imprecision
//This function corrects the lengths
void Cell_two_IS::stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS()
{

    //resizing of the MTOC
    this->MTOC.resize_from_originals(  );


    for( unsigned int micro_index = 0 ; micro_index < this->number_of_microtubules ; micro_index ++ )
    {
        if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 20  )
        {

	    //Checks how many dyneins are on microtubule
	    unsigned int number_of_dynein_1 = this->array_Of_Microtubules[ micro_index ].get_Dynein_points_without_erasing().size();
	    std::vector<Vector3d> motors_for_surface = this->array_Of_Microtubules[ micro_index ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple();
	    //Detached dyneins
	    unsigned int number_of_dynein_2 = 0;
	    if(  this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 20 ) 
	    {
	    	number_of_dynein_2 = this->array_Of_Microtubules[ micro_index ].get_Dynein_points_without_erasing().size();	    
	    }	  
	    //This just check that no dynein was lost during the resizing
	    //It would be a serious error and entire simulation would be stopped   	    
	    if( number_of_dynein_2 + motors_for_surface.size() != number_of_dynein_1 )
	    {
	    	cout<<"number_of_dynein_2 + motors_for_surface.size() != number_of_dynein_1"<<endl;
	    	cout<<"number_of_dynein_1 = "<<number_of_dynein_1<<endl;
	    	cout<<"number_of_dynein_2 = "<<number_of_dynein_2<<endl;
	    	cout<<"motors_for_surface.size() = "<<motors_for_surface.size()<<endl;
	    	cout<<"cell_id = "<<this->cell_id<<endl;
	    	throw("");	    	
	    }
	    //The dyneins are added to the plasma membrane  
	    IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
	    //New points are created
            std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( motors_for_surface.size() , tmp );
            //And added
            this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
        }
        //The same for microtubules attached in the second IS capture-shrinkage
        else if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 40  )
        {
	          unsigned int number_of_dynein_1 = this->array_Of_Microtubules[ micro_index ].get_Dynein_points_without_erasing().size();        
	          std::vector<Vector3d> motors_for_surface = this->array_Of_Microtubules[ micro_index ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple();
		  unsigned int number_of_dynein_2 = 0;
		  if(  this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 40 ) 
		  {
		    	number_of_dynein_2 = this->array_Of_Microtubules[ micro_index ].get_Dynein_points_without_erasing().size();	    
		  }	  
		  if( number_of_dynein_2 + motors_for_surface.size() != number_of_dynein_1 )
		  {
		    	cout<<"number_of_dynein_2 + motors_for_surface.size() != number_of_dynein_1"<<endl;
		    	cout<<"number_of_dynein_1 = "<<number_of_dynein_1<<endl;
		    	cout<<"number_of_dynein_2 = "<<number_of_dynein_2<<endl;
		    	cout<<"motors_for_surface.size() = "<<motors_for_surface.size()<<endl;
		    	cout<<"cell_id = "<<this->cell_id<<endl;
		    	throw("");	    	
		  }	          
	          
	          
	          IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage_2;
                  std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein_2.create_capture_shrinkage_points( motors_for_surface.size() , tmp );
                  this->Capture_Shrinkage_dynein_2.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
        }
        else if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 9  )
        {

            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
            //Stepping of dyneins on the microtubule attached to cortical sliding
            //It returns detached dyneins            
            std::vector<Vector3d> motors_from_micro = this->array_Of_Microtubules[ micro_index ].stepping_detach_real_dynein_abscissa_projection();
            //Contrarily to the case of the capture-shrinkage, both cortical sliding IS are described with one map
            //The reason is that the filaments can pass through one IS and reach the second
	    //The detached dyneins must be returned to the right IS            
            //Some dyneins can be projected outside IS-they have to be returned to the IS

	    //It checks whether the detached points belong to the first IS
	    std::vector<Vector3d> points_inside;
            std::vector<Vector3d> returned_outside_IS = this->IS_Cortic_Sl.control_IS_Cortical_Sliding_points_2_b( motors_from_micro , points_inside );
            //returned_outside_IS definitely do not belong to the IS 
            //points_inside are checked            
	    std::vector<Vector3d> margin_points;
            std::vector<Vector3d> points_inside_inner_ring = this->IS_Cortic_Sl.control_of_inside_and_margins( points_inside , margin_points );
	
	    //Points_inside_inner_ring are projected to the plasma membrane	
	    if( points_inside_inner_ring.size() > 0 )
	    {	
	    	this->project_and_add_points_to_surface( points_inside_inner_ring );
	    }	

	    //margin_points are projected to the plasma membrane in the proximity of the IS, but outside
	    //Such points are returned to the IS 
	    if( margin_points.size() > 0 )
	    {
		ISCorticalSl IS_tmp_1 = this->get_IS_cortical_sliding_first();
		//New points are created	
		std::vector<Vector3d> tmp_points = this->density_surface_dynein.create_cortical_sliding_points( margin_points.size() , IS_tmp_1 );
		//And added
                this->project_and_add_points_to_surface( tmp_points );
            }	



	    //The points that do not belong to the first IS, are checked whether they belong to the second
	    //The procedure is the same like in the case of the first IS
	    std::vector<Vector3d> points_inside_2;
            std::vector<Vector3d> returned_outside_IS_2 = this->IS_Cortic_Sl_2.control_IS_Cortical_Sliding_points_2_b( returned_outside_IS , points_inside_2 );
	    std::vector<Vector3d> margin_points_IS_2;
            std::vector<Vector3d> points_inside_inner_ring_2 = this->IS_Cortic_Sl_2.control_of_inside_and_margins( points_inside_2 , margin_points_IS_2 );
	
	    if( points_inside_inner_ring_2.size() > 0 )
	    {	
	    	this->project_and_add_points_to_surface( points_inside_inner_ring_2 );
	    }	

	    if( margin_points_IS_2.size() > 0 )
	    {
		ISCorticalSl IS_tmp_2 = this->get_IS_cortical_sliding_second();
		std::vector<Vector3d> tmp_points_2 = this->density_surface_dynein.create_cortical_sliding_points( margin_points_IS_2.size() , IS_tmp_2 );
                this->project_and_add_points_to_surface( tmp_points_2 );
            }	


	    //returned_outside_IS_2 do not belong to IS
	    //This is serious imistake
	    //In such a case, error stops the entire simultion
	    if( returned_outside_IS_2.size() > 0 )
	    {	
	    	cout<<"void Cell_two_IS::stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS()"<<endl;
	    	cout<<"returned_outside_IS_2.size() > 0"<<endl;
		throw("");
	    }	

        }
        //Unattached MTs, they are resized
        else if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 0  )
        {
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
        }
        else
        {
            cout<<" this->array_Of_Microtubules[ micro_index ].get_dynein_index() == "<<this->array_Of_Microtubules[ micro_index ].get_dynein_index()<<endl;
            cout<<" void Cell_two_IS::stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS()"<<endl;
	    cout<<"ERROR_ID = 465461664"<<endl;	
	    throw(" ");
        }
    }
}



//This function resizes the cytoskeleton after every step
//The projection operator project forces on hyperplane so that they are perpendicular to the constraints
//Still, the lenghts of the microtubules can change due to numerical imprecision
//This function corrects the lengths
//It is used in the timeDevelopment_two_Cell_calib function - the dyneins are not considered
void Cell_two_IS::stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS_calibration_axis_z()
{
    //resizing of the MTOC
    this->MTOC.resize_from_originals_axis_calibration(  );

    for( unsigned int micro_index = 0 ; micro_index < this->number_of_microtubules ; micro_index ++ )
    {
  
        if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 0  )
        {
            //resizing of microtubules        
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
        }
        else
        {
            cout<<" this->array_Of_Microtubules[ micro_index ].get_dynein_index() == "<<this->array_Of_Microtubules[ micro_index ].get_dynein_index()<<endl;
            cout<<" void Cell_two_IS::stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS_calibration_axis_z()"<<endl;
	    cout<<"ERROR_ID = 311351351312"<<endl;	
	    throw(" ");
        }

    }
}







//In  this function the new dynein attach to the microtubule previously attached in the capture-shrinkage IS
void Cell_two_IS::microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro_two_IS( )
{
    for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++ )
    {
            if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
            {
            	//Dynein attached to MTs in the first IS
                this->microtubule_catch_pair_abscissa_real_dynein_in_IS( microtubule );
            }
            else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 )
            {
            	//Dynein attached to MTs in the second IS            
		this->microtubule_catch_pair_abscissa_real_dynein_in_IS_2( microtubule );
            }
    }

}



void Cell_two_IS::microtubule_catch_pair_abscissa_real_dynein_in_IS_2( unsigned int microtubule )
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


        bool answer = this->IS_Capture_Shrinkage_2.check_IS_capture_segment_control( bead_position , tangent );
        if( ( answer == true ) && ( this->Capture_Shrinkage_dynein_2.get_dynein_points( 1 ).size() > 0 ) )
        {

            std::vector<Vector3d> replace_motors;
            std::vector<Vector3d> motors_to_attach;
            std::vector<Vector3d> dynein_points = this->Capture_Shrinkage_dynein_2.get_dynein_points( 1 );

            for( unsigned int dynein_index = 0 ; dynein_index < dynein_points.size() ; dynein_index ++ )
            {

                Vector3d motor_position = dynein_points[ dynein_index ];
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );

                //double probability = rand_x( 0 , 1 );
		double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                if( probability < Dynein_real::calculate_attachment_probability_per_time_step( distance ) )
                {
		    //cout<<"answer == true "<<endl;
                    motors_to_attach.push_back( motor_position );
                }
                else
                {
                    replace_motors.push_back( motor_position );
                }

            }
		//cout<<"motors_to_attach = "<<motors_to_attach.size()<<endl;
		//cout<<"replace_motors = "<<replace_motors.size()<<endl;

            for( unsigned int counter_id = 0 ; counter_id < motors_to_attach.size() ; counter_id ++ )
            {
		//cout<<"counter_id = "<<endl;
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
		//cout<<"replace_motors"<<endl;
                this->Capture_Shrinkage_dynein_2.erase_vector_points_with_key( 1 );
            }
	    else
	    {
                this->Capture_Shrinkage_dynein_2.set_dynein_points( 1 , replace_motors );
            }


        }

    }


}



//Check, whether unattached MT attaches to capture-shrinkage dynein 
void Cell_two_IS::check_caught_micro_IS_with_real_dynein_2( unsigned int microtubule )
{
     //This is the function for control, whether the microtubule 0 or 9 is caught

    std::uniform_real_distribution<> distribution{ 0 , 1 };
    //This just checks, whether this functions is called on already caught MTs
    if( ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 ) || ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 ) )
    {
	cout<<"this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 ) || ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 "<<endl;
	cout<<" void Cell_two_IS::check_caught_micro_IS_with_real_dynein_2( unsigned int microtubule ) "<<endl;
	unsigned int ERROR_ID = 16516;
	cout<<"ERROR_ID = "<<ERROR_ID<<endl;
	throw("");
    }

    unsigned int dynein_index_to_remember;

    unsigned int lower_index_tmp = 1;
    //Only the last segment can attach - capture-shrinkage acts in end-on fashion
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
	//cout<<"bbbbbbbbbbbbbbbbbbbbbbbbb"<<endl;
        Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d tangent;
        try
        {
            tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );
        }
        catch( int e )
        {
            cout<<"exception"<<endl;
            cout<<"void Cell_two_IS::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )"<<endl;
            throw("");
        }
        //Check whether the segment is in a proximity of the IS        
        bool answer = this->IS_Capture_Shrinkage_2.check_IS_capture_segment_control( bead_position , tangent );
        
        if( answer == true  )
        {
            std::vector<Vector3d> dynein_points = this->Capture_Shrinkage_dynein_2.get_dynein_points( 1 );
            std::vector< Vector3d > replace_motors;


	    //In this cycle the attachment of the points is tested
            for( unsigned int dynein_index = 0 ; dynein_index < dynein_points.size() ; dynein_index ++ )
            {

                Vector3d motor_position = dynein_points[ dynein_index ];
                //Closest point on the segment from the position of the dynein motor                
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                //Distance between the point and the dynein                
                double distance = distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );


		unsigned int number_of_generator = omp_get_thread_num();
    		double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

                if( probability < Dynein_real::calculate_attachment_probability_per_time_step( distance ) )
                {
			dynein_index_to_remember = dynein_index;
			//If  the microtubule is attached to cortical sliding dyneins, they are released 
                        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9  )
                        {
                            std::vector< Vector3d  > points_positions = this->array_Of_Microtubules[ microtubule ].get_Dynein_points_and_erase();
                            this->project_and_add_points_to_surface( points_positions );
                        }



			Vector3d position = motor_position;
			//The position of the motor is chosen as the point, where the microtubule will depolymerize			
                        Vector3d caught_position = this->project_point_on_surface( position );



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
					//cout<<"---40---"<<endl;
					//cout<<distance_new_tangent - this->array_Of_Microtubules[ microtubule ].get_last_Tangent( ).norm()<<endl;
					continue;
				}
				else
				{

				}
			}
         		//Here the depolymeriyzation position is set
                        this->array_Of_Microtubules[ microtubule ].set_IS_position_catching( caught_position  );
                        unsigned int number_of_beads = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1;

                        //substraction of the beads

                        for( unsigned int sub_id = number_of_beads ; sub_id > segment_id ; sub_id -- )
                        {
                        	this->array_Of_Microtubules[ microtubule ].Substract_one_Bead();
                        }
                        this->array_Of_Microtubules[ microtubule ].add_one_final_Bead_catching_position();
			//The growth of the microtubule is stopped	
			this->array_Of_Microtubules[ microtubule ].set_growth_index( 0 );
			this->array_Of_Microtubules[ microtubule ].set_dynein_index( 40 );
			//The length of the microtubule and lenghts are set			
                        this->array_Of_Microtubules[ microtubule ].set_lenght_after_catching( this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule() );
			this->array_Of_Microtubules[ microtubule ].set_lenght_of_tangents("void Cell_two_IS::check_caught_micro_IS_with_real_dynein_2( unsigned int microtubule )");
			break; 
		}
	    }

	}
    }

    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40  )
    {
    	this->array_Of_Microtubules[ microtubule ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple_after_catching();
    }


    //Attaches the point from the previous loop
    //The rest of dyneins are returned to the IS
    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40  )
    {
	
    	std::vector<Vector3d> dynein_points = this->Capture_Shrinkage_dynein_2.get_dynein_points( 1 );
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
   	this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );

 	dynein_points.erase( dynein_points.begin() + dynein_index_to_remember );

    	this->Capture_Shrinkage_dynein_2.set_dynein_points( 1 , dynein_points );

    }


}



//Check whether a new microtubule attaches to capture-shrinkage dynein
void Cell_two_IS::check_if_micro_is_caught_in_IS_two_IS( )
{
    //the condition is here to prevent destruction of cytoskeleton by capture shrinkage 
    for( unsigned int micro_id = 0 ; micro_id < this->get_microtubule_number() ; micro_id ++ )
    {
        if( ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0 ) || ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 9 ) )
        {
            //In the first IS
            this->check_caught_micro_IS_with_real_dynein( micro_id );
        }
        if( ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0 ) || ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 9 ) )
        {
            //In the second IS        
            this->check_caught_micro_IS_with_real_dynein_2( micro_id );
        }
    }

}












//Saving the results 
//Saves results for statistical analysis
void Cell_two_IS::BIG_PRINT_numerical_two_IS( unsigned int time_Step , unsigned int key )
{
   //Saves the center of the MTOC
    this->print_center_MTOC_numerical( time_Step , key );
    //Saves the cortical sliding dynein acting on both IS
    this->print_dynein_cortical_sliding_numerical_two_IS( time_Step , key );
    //Saves the capture-shrinkage dynein acting on both IS
    this->print_dynein_capture_shrinkage_numerical_two_IS( time_Step , key );    
    //Saves the number of microtubules attached in both capture-shrinkage IS 
    this->print_microtubules_capture_shrinkage_numerical_two_IS( time_Step , key );   
    //Saves the number of microtubules shorter than a treshold
    this->print_number_of_shorter_microtubules_two_IS( time_Step , key );  

}


//Saves the number of microtubules shorter than a treshold
void Cell_two_IS::print_number_of_shorter_microtubules_two_IS( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file_1 = path + "NUMBER_OF_SHORTER_MICROTUBULES_" + number + ".txt";

      unsigned int number_of_short_microtubules = this->get_number_of_short_micrtobules();
      ofstream fout_1;
      fout_1.open( name_of_text_file_1.c_str() , std::fstream::app );
      fout_1<<number_of_short_microtubules<<endl;  
      fout_1.close();
}



//Saves the data about unattached and cortical sliding microtubules
//How many microtubules are growing, shrinking, stable and what is the average length
void Cell_two_IS::print_MT_stats( double time_Step , unsigned int key  )
{
	//Microtubules are distinguished by the dynein index determining whether they are attached to dynein
	//Growing index determines their growing, shrinking 
	MatrixXd MT_0 = MatrixXd::Zero( 4 , 1 );
	MatrixXd MT_9 = MatrixXd::Zero( 4 , 1 );

	for( unsigned int index = 0 ; index < this->get_microtubule_number() ;  index ++ )
	{
	        //Unattached microtubules	
		if( this->getMicrotubule( index ).get_dynein_index() == 0 )
		{
			if( this->getMicrotubule( index ).get_growth_index( ) == 0 )
			{
				//The number of stable microtubules is increased			
				MT_0( 0 ) = MT_0( 0 ) + 1;
				//The lenght is added				
				MT_0( 3 ) = MT_0( 3 ) + this->getMicrotubule( index ).get_lenght_of_microtubule_outside_MTOC();
			}
			else if( this->getMicrotubule( index ).get_growth_index( ) == 1 )
			{
				//The number of gwowing microtubules is increased						
				MT_0( 1 ) = MT_0( 1 ) + 1;
				MT_0( 3 ) = MT_0( 3 ) + this->getMicrotubule( index ).get_lenght_of_microtubule_outside_MTOC();
			}
			else if( this->getMicrotubule( index ).get_growth_index( ) == 2 )
			{
				//The number of shrinking microtubules is increased									
				MT_0( 2 ) = MT_0( 2 ) + 1;
				MT_0( 3 ) = MT_0( 3 ) + this->getMicrotubule( index ).get_lenght_of_microtubule_outside_MTOC();
			}
		}
		//The same is done for the microtubule attached to cortical sliding
                if( this->getMicrotubule( index ).get_dynein_index() == 9 )
		{
			if( this->getMicrotubule( index ).get_growth_index( ) == 0 )
			{
				MT_9( 0 ) = MT_9( 0 ) + 1;
				MT_9( 3 ) = MT_9( 3 ) + this->getMicrotubule( index ).get_lenght_of_microtubule_outside_MTOC();
			}
			else if( this->getMicrotubule( index ).get_growth_index( ) == 1 )
			{
				MT_9( 1 ) = MT_9( 1 ) + 1;
				MT_9( 3 ) = MT_9( 3 ) + this->getMicrotubule( index ).get_lenght_of_microtubule_outside_MTOC();
			}
			else if( this->getMicrotubule( index ).get_growth_index( ) == 2 )
			{
				MT_9( 2 ) = MT_9( 2 ) + 1;
				MT_9( 3 ) = MT_9( 3 ) + this->getMicrotubule( index ).get_lenght_of_microtubule_outside_MTOC();
			}
		}


	}
	//The average length is calculated	
	if( ( MT_0( 0 ) + MT_0( 1 ) + MT_0( 2 )  ) >= 1 )
	{
		MT_0( 3 ) = MT_0( 3 ) / ( MT_0( 0 ) + MT_0( 1 ) + MT_0( 2 )  ); 
	}

	if(MT_9( 0 ) + MT_9( 1 ) + MT_9( 2 ) >= 1 )
	{
		MT_9( 3 ) = MT_9( 3 ) / ( MT_9( 0 ) + MT_9( 1 ) + MT_9( 2 )  );
	}


      double time = this->time_clock;
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string( key );


      std::string name_of_text_file = path + "MT_0_" + number + ".txt";
      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<time<<" "<<MT_0( 0 )<< " " << MT_0( 1 ) <<" "<<MT_0( 2 )<<" "<<MT_0( 3 )<<endl;  //write to file
      fout.close();



      std::string name_of_text_file_9 = path + "MT_9_" + number + ".txt";
      fout.open( name_of_text_file_9.c_str() , std::fstream::app );
      fout<<time<<" "<<MT_9( 0 )<< " " << MT_9( 1 ) <<" "<<MT_9( 2 )<<" "<<MT_9( 3 )<<endl;  //write to file
      fout.close();



}	


//Saves the lenght of all microtubules
void Cell_two_IS::print_MT_stats_lenght( double time_Step , unsigned int key  )
{

      std::string path = "./picturesVideos/numerical_results/";
      std::string number_string = std::to_string( key );
      unsigned int time_INT = ( unsigned int ) time_Step;	
      std::string time_string = std::to_string( time_INT );
      std::string name_of_text_file = path + "micro_Lenght_" + time_string + "_"+ number_string + ".txt";


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );

   for( unsigned int microtubule_index = 0 ; microtubule_index < this->number_of_microtubules ; microtubule_index ++ )
   {
      fout<< this->getMicrotubule( microtubule_index ).get_dynein_index() + this->getMicrotubule( microtubule_index ).get_growth_index( )<<" "<< this->getMicrotubule( microtubule_index ).get_lenght_of_microtubule_outside_MTOC()<<endl;
   }
   fout.close();




}	





//Saves the number of microtubules attached in both capture-shrinkage IS 
void Cell_two_IS::print_microtubules_capture_shrinkage_numerical_two_IS( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file_1 = path + "IS_CAPTURE_SHRINKAGE_MICROTUBULES_1_" + number + ".txt";
      std::string name_of_text_file_2 = path + "IS_CAPTURE_SHRINKAGE_MICROTUBULES_2_" + number + ".txt";

      unsigned int number_of_dynein_1 = this->get_number_of_microtubules_micro_capture_20();
      unsigned int number_of_dynein_2 = this->get_number_of_microtubules_micro_capture_40();

      ofstream fout_1;
      fout_1.open( name_of_text_file_1.c_str() , std::fstream::app );
      fout_1<<number_of_dynein_1<<endl;  //write to file
      fout_1.close();

      ofstream fout_2;
      fout_2.open( name_of_text_file_2.c_str() , std::fstream::app );
      fout_2<<number_of_dynein_2<<endl;  //write to file
      fout_2.close();
}


//Saves the number of capture-shrinkage dyneins from both IS acting on microtubules
void Cell_two_IS::print_dynein_capture_shrinkage_numerical_two_IS( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file_1 = path + "IS_CAPTURE_SHRINKAGE_1_" + number + ".txt";
      std::string name_of_text_file_2 = path + "IS_CAPTURE_SHRINKAGE_2_" + number + ".txt";

      unsigned int number_of_dynein_1 = this->get_number_of_dynein_motors_capture_shrinkage_1();
      unsigned int number_of_dynein_2 = this->get_number_of_dynein_motors_capture_shrinkage_2();

      ofstream fout_1;
      fout_1.open( name_of_text_file_1.c_str() , std::fstream::app );
      fout_1<<number_of_dynein_1<<endl;  //write to file
      fout_1.close();

      ofstream fout_2;
      fout_2.open( name_of_text_file_2.c_str() , std::fstream::app );
      fout_2<<number_of_dynein_2<<endl;  //write to file
      fout_2.close();
}



//This functions saves the number of cortical sliding dyneins from both IS acting on microtubules 
void Cell_two_IS::print_dynein_cortical_sliding_numerical_two_IS( double time_Step , unsigned int key )
{
      //The number of cortical sliding dynein is saved
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "IS_CORTICAL_SLIDING_" + number + ".txt";
      unsigned int number_of_dynein = this->get_number_of_dynein_motors_cortical_sliding_micro();


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<number_of_dynein<<endl;  //write to file
      fout.close();

      std::vector<Vector3d> dynein_points = this->density_surface_dynein.get_all_dynein_point();

      //The cortical sliding dyneins from IS1
      unsigned int counter_IS_1 = 0;
      //The cortical sliding dyneins from IS2
      unsigned int counter_IS_2 = 0;

     for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++ )
     {
	   if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
	   {
	   	std::vector<Vector3d> dynein_points = this->array_Of_Microtubules[ microtubule ].get_dynein_points_Cortical_Sliding();
	   	//Dyneins from both IS are counted	   	
           	for( unsigned int index = 0 ; index < dynein_points.size() ; index ++ )
           	{
			Vector3d dynein_position = dynein_points[ index ];
			bool answer_1 = this->IS_Cortic_Sl.control_IS_Cortical_Sliding_points( dynein_position );
			bool answer_2 = this->IS_Cortic_Sl_2.control_IS_Cortical_Sliding_points( dynein_position );
			if ( answer_1 == true )
			{
				counter_IS_1++;
			}
			else if ( answer_2 == true )
			{
				counter_IS_2++;
			}
           	}
	   }
     }

      std::string name_of_text_file_1 = path + "IS_CORTICAL_SLIDING_1_" + number + ".txt";
      std::string name_of_text_file_2 = path + "IS_CORTICAL_SLIDING_2_" + number + ".txt";

      ofstream fout_1;
      fout_1.open( name_of_text_file_1.c_str() , std::fstream::app );
      fout_1<<counter_IS_1<<endl;  //write to file
      fout_1.close();

      ofstream fout_2;
      fout_2.open( name_of_text_file_2.c_str() , std::fstream::app );
      fout_2<<counter_IS_2<<endl;  //write to file
      fout_2.close();
}



//Returns the number of capture-shrinkage dyneins in the IS 1 attached to the microtubules
unsigned int Cell_two_IS::get_number_of_dynein_motors_capture_shrinkage_1()
{

    unsigned int sum_of_motor_in_micro = 0;
    for( unsigned int motor_id = 0 ; motor_id < this->number_of_microtubules; motor_id ++ )
    {
        if( this->array_Of_Microtubules[ motor_id ].get_dynein_index() == 20 )
        {
            unsigned int number_of_dynein_in_micro = this->array_Of_Microtubules[ motor_id ].get_Dynein_abscissa().size();
            sum_of_motor_in_micro = sum_of_motor_in_micro + number_of_dynein_in_micro;
        }
    }
    return sum_of_motor_in_micro;
}



//Returns the number of capture-shrinkage dyneins in the IS 2 attached to the microtubules   
unsigned int Cell_two_IS::get_number_of_dynein_motors_capture_shrinkage_2()
{

    unsigned int sum_of_motor_in_micro = 0;
    for( unsigned int motor_id = 0 ; motor_id < this->number_of_microtubules; motor_id ++ )
    {
        if( this->array_Of_Microtubules[ motor_id ].get_dynein_index() == 40 )
        {
            unsigned int number_of_dynein_in_micro = this->array_Of_Microtubules[ motor_id ].get_Dynein_abscissa().size();
            sum_of_motor_in_micro = sum_of_motor_in_micro + number_of_dynein_in_micro;
        }
    }
    return sum_of_motor_in_micro;
}



//Returns the number of microtubules attached to capture-shrinkage dyneins in IS1 
unsigned int Cell_two_IS::get_number_of_microtubules_micro_capture_20()
{
	unsigned int counter_20 = 0;
	for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules; microtubule ++  )
	{
	        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
	        {
	            counter_20 = counter_20++;
	        }

	}
	return counter_20;
}

//Returns the number of microtubules attached to capture-shrinkage dyneins in IS2
unsigned int Cell_two_IS::get_number_of_microtubules_micro_capture_40()
{
	unsigned int counter_40 = 0;
	for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules; microtubule ++  )
	{
	        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 )
	        {
	            counter_40 = counter_40++;
	        }

	}
	return counter_40;
}




//This function gives the number of microtubules, whose number of beads is smaller than a treshold
//The purpose of the function is only a control
unsigned int Cell_two_IS::get_number_of_short_micrtobules()
{
	unsigned int counter_shorter = 0;
	for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules; microtubule ++  )
	{
	        if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() < ( sim_of_Cell::MicrotubulePoints - sim_of_Cell::REDUCTION - 2 )  )
	        {
	            counter_shorter++;
	        }

	}
	return counter_shorter;
}



//This functions control, whether microtubules with no dynein were properly detached from the IS
void Cell_two_IS::control_dynein_detachment()
{
	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		//kontroluju, jestli tam jsou dyneiny
		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 40 )
		{
			//If the microtubule has no dynein, it detaches 
			if( this->array_Of_Microtubules[ micro_id ].get_number_of_dynein_points_IS() == 0 )
			{				
				IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage_2;
 				std::vector<Vector3d> vectors_to_return = this->array_Of_Microtubules[ micro_id ].detach_from_IS_2();
				//This sets the length of microtubule after the detachment
				if( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() <  this->array_Of_Microtubules[ micro_id ].get_lenght_of_micro_0_outside_MTOC(  ) )
				{
					this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );
				}
				else
				{
					//Some MTs can prolong  due to the attachment to dynein, as the depolymerization point is set at the position of the dynein
					//In such a case, their length is controlled and possibly set back 				
					if( this->array_Of_Microtubules[ micro_id ].control_of_captured_MT_lenght_against_prolongation() == true )
					{
						this->array_Of_Microtubules[ micro_id ].shorten_MT_due_to_numerical_prolongation_during_attachment( );
					}

				}

				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0();
                                this->array_Of_Microtubules[ micro_id ].set_growth_index( 0 );
                                //Simple control, whether there were really no dyneins
                                //The opposite case would be a mistake and the simulation would be stopped
				if( vectors_to_return.size() > 0 )
				{
					cout<<"vectors_to_return.size() > 0"<<endl;
					cout<<"ERROR_ID = 5541413113"<<endl;
					cout<<"void Cell_two_IS::control_dynein_detachment()"<<endl;
					throw("");
				}

			}
		}
		//The same is done for the microtubule attached to the first IS		
		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 20 )
		{

			if( this->array_Of_Microtubules[ micro_id ].get_number_of_dynein_points_IS() == 0 )
			{
				IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
 				std::vector<Vector3d> vectors_to_return = this->array_Of_Microtubules[ micro_id ].detach_from_IS_2();
				if( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() <  this->array_Of_Microtubules[ micro_id ].get_lenght_of_micro_0_outside_MTOC(  ) )
				{
					this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );
				}
				else
				{

					if( this->array_Of_Microtubules[ micro_id ].control_of_captured_MT_lenght_against_prolongation() == true )
					{
						this->array_Of_Microtubules[ micro_id ].shorten_MT_due_to_numerical_prolongation_during_attachment( );
					}

				}

				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_growth_index( 0 );
				if( vectors_to_return.size() > 0 )
				{
					cout<<"vectors_to_return.size() > 0"<<endl;
					cout<<"ERROR_ID = 55414131444646613"<<endl;
					cout<<"void Cell_two_IS::control_dynein_detachment()"<<endl;
					throw("");
				}


			}
		}
	}
}





//Controls the prolongation of microtubules
void Cell_two_IS::control_lenght_prolongation_capture()
{
	//Currently unused
	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		//kontroluju, jestli tam jsou dyneiny
		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 20 )
		{
			if( this->array_Of_Microtubules[ micro_id ].control_of_captured_MT_lenght_against_prolongation() == true )
			{
 				std::vector<Vector3d> vectors_to_return = this->array_Of_Microtubules[ micro_id ].detach_from_IS_2();
				this->array_Of_Microtubules[ micro_id ].shorten_MT_due_to_numerical_prolongation_during_attachment( );
				this->array_Of_Microtubules[ micro_id ].set_growth_index( 0 );
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::control_lenght_prolongation_capture()");	
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC(  ) );




				IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
            			std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( vectors_to_return.size() , tmp );
            			this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );

			}
		}
		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 40 )
		{
			if( this->array_Of_Microtubules[ micro_id ].control_of_captured_MT_lenght_against_prolongation() == true )
			{
 				std::vector<Vector3d> vectors_to_return = this->array_Of_Microtubules[ micro_id ].detach_from_IS_2();
				this->array_Of_Microtubules[ micro_id ].shorten_MT_due_to_numerical_prolongation_during_attachment( );
				this->array_Of_Microtubules[ micro_id ].set_growth_index( 0 );
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::control_lenght_prolongation_capture()");	
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC(  ) );



				IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage_2;
            			std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein_2.create_capture_shrinkage_points( vectors_to_return.size() , tmp );
            			this->Capture_Shrinkage_dynein_2.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );

			}
		}
	}
}


















void Cell_two_IS::all_micro_growth()
{
	std::uniform_real_distribution<> distribution{ 0 , 1 };
	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		if( this->array_Of_Microtubules[ micro_id ].getNumberOfPoints() >= sim_of_Cell::MicrotubulePoints - sim_of_Cell::REDUCTION )
		{
			continue;
		}

		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0  )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index( ) == 0 )
			{
		    		unsigned int number_of_generator = omp_get_thread_num();
		    		double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				double probability_rescue_per_time_step = dynamic_instability::rescue_rate * sim_of_Cell::time_Step;
				if( probability <= probability_rescue_per_time_step )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
			}
		}
	}
	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0  )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index( ) == 1 )
			{

				this->array_Of_Microtubules[ micro_id ].microtubule_growth_2();
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::all_micro_growth()");
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	
				if( this->array_Of_Microtubules[ micro_id ].getNumberOfPoints() == sim_of_Cell::MicrotubulePoints )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 0 );
				}
			}

		}
	}



}



void Cell_two_IS::set_dynamic_instability()
{
	//get_dynein_index - 0 default state
	//                   1 growing
	//		     2 shrinking	

	unsigned int number_of_generator = omp_get_thread_num();
	double ratio_of_rates = dynamic_instability::rescue_rate / ( dynamic_instability::catastrophe_rate + dynamic_instability::rescue_rate );



	std::uniform_real_distribution<> distribution{ 0 , 1 };
	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 0 )
		{
			cout<<" this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 0 "<<endl;
			cout<<"void Cell_two_IS::set_dynamic_instability()"<<endl;
			cout<<"ERROR_ID = 315613156"<<endl;
			throw("");			
		}

		if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 0  )
		{
		    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
			if( probability < ratio_of_rates )
			{
				this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
			}
			else 
			{
				this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
			}
		}
	

	}

}

















/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




MatrixXd Cell_two_IS::dynamic_instability_MT_grow_linear(  )
{


	unsigned int Number_of_MT = this->get_microtubule_number();
	//unsigned int Number_of_MT = 1e5;
	MatrixXd Matrix_MTs = MatrixXd::Zero( Number_of_MT , 2 );	

	//initialization of length
	double initial_lenght = sim_of_Cell::PI * this->a_axis;//sim_of_Cell::PI * this->a_axis
	for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
	{
		Matrix_MTs( micro_index , 0 ) = initial_lenght;
		Matrix_MTs( micro_index , 1 ) = 1.0;	
	}

        std::uniform_real_distribution<> distribution{ 0 , 1 };
        unsigned int number_of_generator = omp_get_thread_num();


////////////////////////////////////////////////////////////////////

	//dynamic instability process
	//unsigned int number_of_time_steps = (unsigned int)( dynamic_instability::Time_Of_initial_process / sim_of_Cell::time_Step );	
	unsigned int number_of_time_steps = (unsigned int)( dynamic_instability::Time_Of_initial_process / dynamic_instability::Time_Step_Of_in_pr );
	for( unsigned int time_step_id = 0 ; time_step_id <  number_of_time_steps ; time_step_id ++ )
	{

		//grow and shrinkage
		for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
		{
			if( abs( Matrix_MTs( micro_index , 1 ) - 1.0 ) < 1e-6 )
			{
				Matrix_MTs( micro_index , 0 ) = Matrix_MTs( micro_index , 0 ) + dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::polymerization_constant;
			}
			else if( abs( Matrix_MTs( micro_index , 1 ) - 2.0 ) < 1e-6 )
			{
				Matrix_MTs( micro_index , 0 ) = Matrix_MTs( micro_index , 0 ) - dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::shrinking_constant;				
			}
			else
			{
				cout<<"void Cell_two_IS::dynamic_instability_MT_grow(  )"<<endl;
				cout<<"Matrix_MTs( micro_index , 1 ) =  "<<Matrix_MTs( micro_index , 1 )<<endl;
				cout<<"ERROR_ID = 653131313513"<<endl;
				throw("");
			}
			if( Matrix_MTs( micro_index , 0 )  <= dynamic_instability::coefficient_4  )
			{
				Matrix_MTs( micro_index , 1 ) = 1.0;
			}
			if( Matrix_MTs( micro_index , 0 )  >= dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  )
			{
				Matrix_MTs( micro_index , 1 ) = 2.0;
			}


		}		


		//growing and shrinking
		for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ;  micro_index ++ )
		{
			if( abs( Matrix_MTs( micro_index , 1 ) - 1.0 ) < 1e-6 )
			{
				double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability <  dynamic_instability::catastrophe_rate * dynamic_instability::Time_Step_Of_in_pr )
				{
					Matrix_MTs( micro_index , 1 ) = 2;
				}
			}
			else if( abs( Matrix_MTs( micro_index , 1 ) - 2.0 ) < 1e-6 )
			{
				double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability <  dynamic_instability::rescue_rate * dynamic_instability::Time_Step_Of_in_pr )
				{
					Matrix_MTs( micro_index , 1 ) = 1;
				}
			}
		}	
	
	}
	//controling for

	for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
	{
		double length = Matrix_MTs( micro_index , 0 );
		if( length < dynamic_instability::coefficient_4 - 2.0 * dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::shrinking_constant )
		{
			cout<<" length < dynamic_instability::coefficient_4 - 2.0 * dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::shrinking_constant  "<<endl;
			cout<<"MatrixXd Cell_two_IS::dynamic_instability_MT_grow_linear(  )"<<endl;
			cout<<"ERROR_ID = 67684646341351361313"<<endl;
			throw("");	
		} 
		else if( length > dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  + 2.0 * dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::polymerization_constant  )
		{
			cout<<" length > dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  + 2.0 * dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  "<<endl;
			cout<<"MatrixXd Cell_two_IS::dynamic_instability_MT_grow_linear(  )"<<endl;
			cout<<"ERROR_ID = 5456133125143125431351"<<endl;
			throw("");
		}
	}
	//narovnani
	for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
	{
		double length = Matrix_MTs( micro_index , 0 );
		if( length < dynamic_instability::coefficient_4 )
		{
			cout<<"----------------------------------"<<endl;
			cout<<"dynamic_instability::coefficient_4 = "<<dynamic_instability::coefficient_4<<endl;
			Matrix_MTs( micro_index , 0 ) = dynamic_instability::coefficient_4 + dynamic_instability::polymerization_constant * 0.5; //increase for half a second
		} 
		else if( length > dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  )
		{
			Matrix_MTs( micro_index , 0 ) = dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI - 0.5 * dynamic_instability::shrinking_constant;//decrease for half a second
		}
	}

	for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
	{
		double length = Matrix_MTs( micro_index , 0 );
		if( length < dynamic_instability::coefficient_4 )
		{
			cout<<" length < dynamic_instability::coefficient_4 - 2.0 * dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::shrinking_constant  "<<endl;
			cout<<"MatrixXd Cell_two_IS::dynamic_instability_MT_grow_linear(  )"<<endl;
			cout<<"ERROR_ID = 67684646341351361313"<<endl;
			throw("");	
		} 
		else if( length > dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  )
		{
			cout<<" length > dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  + 2.0 * dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  "<<endl;
			cout<<"MatrixXd Cell_two_IS::dynamic_instability_MT_grow_linear(  )"<<endl;
			cout<<"ERROR_ID = 5456133125143125431351"<<endl;
			throw("");
		}
	}

         
        char name_of_text_file [100];
	sprintf ( name_of_text_file , "./picturesVideos/micro_Lenght_linear_distribution.txt" );

	FILE *out;
	out = fopen( name_of_text_file , "w");

	for( unsigned int row = 0 ; row < Matrix_MTs.rows() ; row ++ )
	{
		//cout<<"row = "<<row<<endl;
		fprintf(out,"%20.10f %10.5f\n",Matrix_MTs( row , 0 )  , Matrix_MTs( row ,  1 ) );
	}
        fclose( out );
	return Matrix_MTs;
}







void Cell_two_IS::dynamic_instability_linear()
{
	//get_dynein_index - 0 default state
	//                   1 growing
	//		     2 shrinking	


	double probability_rescue_per_time_step = dynamic_instability::rescue_rate * sim_of_Cell::time_Step;
	double probability_catastrophe_per_time_step = dynamic_instability::catastrophe_rate * sim_of_Cell::time_Step;
	unsigned int number_of_generator = omp_get_thread_num();



	std::uniform_real_distribution<> distribution{ 0 , 1 };
	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		if( ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 0 ) && ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 9 ) )
		{
			continue;
		}

		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0 )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 0  )
			{

			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

				if( probability < probability_rescue_per_time_step )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
				else if( ( probability >= probability_rescue_per_time_step ) && ( probability < probability_rescue_per_time_step + probability_catastrophe_per_time_step )  )
				{

					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
				else
				{

				}

			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{
			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_catastrophe_per_time_step )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
				if( this->array_Of_Microtubules[ micro_id ].get_lenght_of_micro_0_outside_MTOC() > dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{
			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_rescue_per_time_step )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
				if( this->array_Of_Microtubules[ micro_id ].get_lenght_of_micro_0_outside_MTOC() < dynamic_instability::coefficient_4 )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}


			}
		}
		else if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 9  )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 0  )
			{

			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

				if( probability < probability_rescue_per_time_step )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
				else if( ( probability >= probability_rescue_per_time_step ) && ( probability < probability_rescue_per_time_step + probability_catastrophe_per_time_step )  )
				{

					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
				else
				{

				}

			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{
			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_catastrophe_per_time_step )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
				if( this->array_Of_Microtubules[ micro_id ].get_lenght_of_micro_0_outside_MTOC() > dynamic_instability::coefficient_3 * this->a_axis * sim_of_Cell::PI  )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{
			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_rescue_per_time_step )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
				if( this->array_Of_Microtubules[ micro_id ].get_lenght_of_micro_0_outside_MTOC() < dynamic_instability::coefficient_4 )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}


			}
		}

		//this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
	}




	//RUST---------------------------------------------------------------------------------


	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		if( ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 0 ) && ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 9 ) )
		{
			continue;
		}

		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0  )
		{

			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{

				this->array_Of_Microtubules[ micro_id ].microtubule_growth_MT_0_9();
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	

			}

			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{
				this->array_Of_Microtubules[ micro_id ].microtubule_shrinkage_MT_0_9();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	

			}

		}






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 9  )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{
				this->array_Of_Microtubules[ micro_id ].microtubule_growth_MT_0_9();
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	
			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{
				    //cout<< 9 << "this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2  "<<endl;
				    this->array_Of_Microtubules[ micro_id ].microtubule_shrinkage_MT_0_9();
				    this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				    this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9();
				    std::vector< Vector3d > motors_from_micro = this->array_Of_Microtubules[ micro_id ].detachmet_due_to_shrinkage_of_MT_9( );
				    this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				    this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	

				    if( motors_from_micro.size() == 0 )
                  		    {
					continue;
				    }	
	
				    std::vector<Vector3d> points_inside;
				    std::vector<Vector3d> returned_outside_IS = this->IS_Cortic_Sl.control_IS_Cortical_Sliding_points_2_b( motors_from_micro , points_inside );
				    std::vector<Vector3d> margin_points;
				    std::vector<Vector3d> points_inside_inner_ring = this->IS_Cortic_Sl.control_of_inside_and_margins( points_inside , margin_points );
				
				    if( points_inside_inner_ring.size() > 0 )
				    {	
					//cout<<"points_inside_inner_ring.size() = "<<points_inside_inner_ring.size()<<endl;
				    	this->project_and_add_points_to_surface( points_inside_inner_ring );
				    }	

				    if( margin_points.size() > 0 )
				    {
					//cout<<"margin_points.size() = "<<margin_points.size()<<endl;
					//zlabava
					ISCorticalSl IS_tmp_1 = this->get_IS_cortical_sliding_first();
					std::vector<Vector3d> margin_points_2 = this->density_surface_dynein.create_cortical_sliding_points( margin_points.size() , IS_tmp_1 );
					//std::vector<Vector3d> margin_points_2 = this->IS_Cortic_Sl.create_IS_Cortical_Sliding_points( margin_points.size() );
					this->project_and_add_points_to_surface( margin_points_2 );
				    }	




				    std::vector<Vector3d> points_inside_2;
				    std::vector<Vector3d> returned_outside_IS_2 = this->IS_Cortic_Sl_2.control_IS_Cortical_Sliding_points_2_b( returned_outside_IS , points_inside_2 );
				    std::vector<Vector3d> margin_points_IS_2;
				    std::vector<Vector3d> points_inside_inner_ring_2 = this->IS_Cortic_Sl_2.control_of_inside_and_margins( points_inside_2 , margin_points_IS_2 );
				
				    if( points_inside_inner_ring_2.size() > 0 )
				    {	
				    	this->project_and_add_points_to_surface( points_inside_inner_ring_2 );
				    }	

				    if( margin_points_IS_2.size() > 0 )
				    {
					ISCorticalSl IS_tmp_2 = this->get_IS_cortical_sliding_second();      
					std::vector<Vector3d> margin_points_2 = this->density_surface_dynein.create_cortical_sliding_points( margin_points_IS_2.size() , IS_tmp_2 );
					//std::vector<Vector3d> margin_points_IS_2_2 = this->IS_Cortic_Sl_2.create_IS_Cortical_Sliding_points( margin_points_IS_2.size() );
					this->project_and_add_points_to_surface( margin_points_2 );
				    }	

				    if( returned_outside_IS_2.size() > 0 )
				    {	
				    	cout<<"void Cell_two_IS::dynamic_instability()"<<endl;
				    	cout<<"returned_outside_IS_2.size() > 0"<<endl;
					throw("");
				    }	



				    cout<< 9 << "-----------------------------------------  "<<endl;
			
			}
		}



	}

}











MatrixXd Cell_two_IS::dynamic_instability_MT_grow(  )
{


	unsigned int Number_of_MT = this->get_microtubule_number();
	Number_of_MT = 5e5;
	MatrixXd Matrix_MTs = MatrixXd::Zero( Number_of_MT , 2 );	

	//initialization of length
	double initial_lenght = sim_of_Cell::PI * this->a_axis;//sim_of_Cell::PI * this->a_axis
	for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
	{
		Matrix_MTs( micro_index , 0 ) = initial_lenght;
		Matrix_MTs( micro_index , 1 ) = 1.0;	
		//cout<<Matrix_MTs( micro_index , 1 )<<endl;
	}



        std::uniform_real_distribution<> distribution{ 0 , 1 };
        unsigned int number_of_generator = omp_get_thread_num();


	//dynamic instability process
	//unsigned int number_of_time_steps = (unsigned int)( dynamic_instability::Time_Of_initial_process / sim_of_Cell::time_Step );	
	unsigned int number_of_time_steps = (unsigned int)( dynamic_instability::Time_Of_initial_process / dynamic_instability::Time_Step_Of_in_pr );
	for( unsigned int time_step_id = 0 ; time_step_id <  number_of_time_steps ; time_step_id ++ )
	{

		//grow and shrinkage
		for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
		{
			if( abs( Matrix_MTs( micro_index , 1 ) - 1.0 ) < 1e-6 )
			{
				Matrix_MTs( micro_index , 0 ) = Matrix_MTs( micro_index , 0 ) + dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::polymerization_constant;
			}
			else if( abs( Matrix_MTs( micro_index , 1 ) - 2.0 ) < 1e-6 )
			{
				Matrix_MTs( micro_index , 0 ) = Matrix_MTs( micro_index , 0 ) - dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::shrinking_constant;				
			}
			else
			{
				cout<<"void Cell_two_IS::dynamic_instability_MT_grow(  )"<<endl;
				cout<<"Matrix_MTs( micro_index , 1 ) =  "<<Matrix_MTs( micro_index , 1 )<<endl;
				cout<<"ERROR_ID = 653131313513"<<endl;
				throw("");
			}
			if( Matrix_MTs( micro_index , 0 )  <=  dynamic_instability::coefficient_4 )
			{
				Matrix_MTs( micro_index , 1 ) = 1.0;
			}


		}		


		//growing and shrinking
		for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ;  micro_index ++ )
		{
			if( abs( Matrix_MTs( micro_index , 1 ) - 1.0 ) < 1e-6 )
			{
				double Lenght_tmp = Matrix_MTs( micro_index , 0 );
				double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability <  dynamic_instability::rate_based_on_L( Lenght_tmp ) * dynamic_instability::Time_Step_Of_in_pr )
				{
					Matrix_MTs( micro_index , 1 ) = 2;
				}
			}
			else if( abs( Matrix_MTs( micro_index , 1 ) - 2.0 ) < 1e-6 )
			{
				double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability <  dynamic_instability::rescue_rate * dynamic_instability::Time_Step_Of_in_pr )
				{
					Matrix_MTs( micro_index , 1 ) = 1;
				}
			}
		}	
	
	}


         
        char name_of_text_file [100];
	sprintf ( name_of_text_file , "./picturesVideos/numerical_results/micro_Lenght_distribution.txt" );

	FILE *out;
	out = fopen( name_of_text_file , "w");

	for( unsigned int row = 0 ; row < Matrix_MTs.rows() ; row ++ )
	{
		//cout<<"row = "<<row<<endl;
		fprintf(out,"%20.10f %10.5f\n",Matrix_MTs( row , 0 )  , Matrix_MTs( row ,  1 ) );
	}
        fclose( out );
	return Matrix_MTs;
}






void Cell_two_IS::dynamic_instability()
{
	//get_dynein_index - 0 default state
	//                   1 growing
	//		     2 shrinking	


	double probability_rescue_per_time_step = dynamic_instability::rescue_rate * sim_of_Cell::time_Step;
	double probability_catastrophe_per_time_step = dynamic_instability::catastrophe_rate * sim_of_Cell::time_Step;
	unsigned int number_of_generator = omp_get_thread_num();



	std::uniform_real_distribution<> distribution{ 0 , 1 };
	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		if( ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 0 ) && ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 9 ) )
		{
			continue;
		}

		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0  )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 0  )
			{

			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

				double L_outside_of_MTOC = this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC();
				double probability_catastrophe_depending_on_lenght = dynamic_instability::rate_based_on_L( L_outside_of_MTOC );
				double probability_catastrophe_depending_on_lenght_per_time_step = probability_catastrophe_depending_on_lenght * sim_of_Cell::time_Step;

				if( probability < probability_rescue_per_time_step )
				{

					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
				else if( ( probability >= probability_rescue_per_time_step ) && ( probability < probability_rescue_per_time_step + probability_catastrophe_depending_on_lenght_per_time_step )  )
				{

					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
				else
				{

				}

			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{
				//tady to menim v zavislosti na delce
				/*
			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_catastrophe_per_time_step )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
				*/
				//tady to menim v zavislosti na delce
				
				double L_outside_of_MTOC = this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC();
				double probability_catastrophe_depending_on_lenght = dynamic_instability::rate_based_on_L( L_outside_of_MTOC );
				double probability_catastrophe_depending_on_lenght_per_time_step = probability_catastrophe_depending_on_lenght * sim_of_Cell::time_Step;

			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_catastrophe_depending_on_lenght_per_time_step )
				{
					//cout<<"LLLLLLLLLLLLLLLLLLLLLL"<<endl;
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{
			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_rescue_per_time_step )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
				if( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() < dynamic_instability::coefficient_4  )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}


			}
		}
		else if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 9  )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 0  )
			{

			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_rescue_per_time_step )
				{
					//cout<<"LLLLLLLLLLLLLLLLLLLLLL"<<endl;
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
				else if( ( probability >= probability_rescue_per_time_step ) && ( probability < probability_rescue_per_time_step + probability_catastrophe_per_time_step )  )
				{
					//cout<<"DDDDDDDDDDDDDDDDDDDDD"<<endl;
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
				else
				{

				}

			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{
				double L_outside_of_MTOC = this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC();
				double probability_catastrophe_depending_on_lenght = dynamic_instability::rate_based_on_L( L_outside_of_MTOC );
				double probability_catastrophe_depending_on_lenght_per_time_step = probability_catastrophe_depending_on_lenght * sim_of_Cell::time_Step;

			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_catastrophe_depending_on_lenght_per_time_step )
				{
					//cout<<"9  LLLLLLLLLLLLLLLLLLLLLL"<<endl;
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
				}
			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{
			    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability < probability_rescue_per_time_step )
				{
					//cout<<"qwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww"<<endl;
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
				if( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() < dynamic_instability::coefficient_4  )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
				}
			}
		}


	}




	//RUST---------------------------------------------------------------------------------











	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		if( ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 0 ) && ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 9 ) )
		{
			continue;
		}

		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0  )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{
				//cout<< 0 << "this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1  "<<endl;
				this->array_Of_Microtubules[ micro_id ].microtubule_growth_MT_0_9();
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability()");
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	
			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{
				//cout<< 0 << "this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2  "<<endl;
				this->array_Of_Microtubules[ micro_id ].microtubule_shrinkage_MT_0_9();
				if( ( this->array_Of_Microtubules[ micro_id ].getNumberOfPoints() == 3 ) && ( this->array_Of_Microtubules[ micro_id ].get_last_Tangent( ).norm() < sim_of_Cell::resting_distance / 2.0 ) )
				{
					this->array_Of_Microtubules[ micro_id ].set_growth_index( 0 );
				}


				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability()");
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	
			}

		}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 9  )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{
				//cout<< 9 << "this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1  "<<endl;
				this->array_Of_Microtubules[ micro_id ].microtubule_growth_MT_0_9();
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability()");
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	
			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{
				    //cout<< 9 << "this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2  "<<endl;
				    this->array_Of_Microtubules[ micro_id ].microtubule_shrinkage_MT_0_9();
				    this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9();
				    std::vector< Vector3d > motors_from_micro = this->array_Of_Microtubules[ micro_id ].detachmet_due_to_shrinkage_of_MT_9( );
				    this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability()");
				    this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	

				    if( motors_from_micro.size() == 0 )
                  		    {
					continue;
				    }	
	
				    std::vector<Vector3d> points_inside;
				    std::vector<Vector3d> returned_outside_IS = this->IS_Cortic_Sl.control_IS_Cortical_Sliding_points_2_b( motors_from_micro , points_inside );
				    std::vector<Vector3d> margin_points;
				    std::vector<Vector3d> points_inside_inner_ring = this->IS_Cortic_Sl.control_of_inside_and_margins( points_inside , margin_points );
				
				    if( points_inside_inner_ring.size() > 0 )
				    {	
					//cout<<"points_inside_inner_ring.size() = "<<points_inside_inner_ring.size()<<endl;
				    	this->project_and_add_points_to_surface( points_inside_inner_ring );
				    }	

				    if( margin_points.size() > 0 )
				    {
					//cout<<"margin_points.size() = "<<margin_points.size()<<endl;
					//zlabava
					ISCorticalSl IS_tmp_1 = this->get_IS_cortical_sliding_first();
					std::vector<Vector3d> margin_points_2 = this->density_surface_dynein.create_cortical_sliding_points( margin_points.size() , IS_tmp_1 );
					//std::vector<Vector3d> margin_points_2 = this->IS_Cortic_Sl.create_IS_Cortical_Sliding_points( margin_points.size() );
					this->project_and_add_points_to_surface( margin_points_2 );
				    }	




				    std::vector<Vector3d> points_inside_2;
				    std::vector<Vector3d> returned_outside_IS_2 = this->IS_Cortic_Sl_2.control_IS_Cortical_Sliding_points_2_b( returned_outside_IS , points_inside_2 );
				    std::vector<Vector3d> margin_points_IS_2;
				    std::vector<Vector3d> points_inside_inner_ring_2 = this->IS_Cortic_Sl_2.control_of_inside_and_margins( points_inside_2 , margin_points_IS_2 );
				
				    if( points_inside_inner_ring_2.size() > 0 )
				    {	
				    	this->project_and_add_points_to_surface( points_inside_inner_ring_2 );
				    }	

				    if( margin_points_IS_2.size() > 0 )
				    {
					ISCorticalSl IS_tmp_2 = this->get_IS_cortical_sliding_second();      
					std::vector<Vector3d> margin_points_2 = this->density_surface_dynein.create_cortical_sliding_points( margin_points_IS_2.size() , IS_tmp_2 );
					//std::vector<Vector3d> margin_points_IS_2_2 = this->IS_Cortic_Sl_2.create_IS_Cortical_Sliding_points( margin_points_IS_2.size() );
					this->project_and_add_points_to_surface( margin_points_2 );
				    }	

				    if( returned_outside_IS_2.size() > 0 )
				    {	
				    	cout<<"void Cell_two_IS::dynamic_instability()"<<endl;
				    	cout<<"returned_outside_IS_2.size() > 0"<<endl;
					throw("");
				    }	



				    cout<< 9 << "-----------------------------------------  "<<endl;
			
			}
		}



	}

}



//This function has two goals: get the probability distribution of microtubule lengths
//It also used to simulate the effects of dynamic instability, so that their length correspond to the probability distribution - length adjustments before the beginning of the simulation
MatrixXd Cell_two_IS::dynamic_instability_MT_grow_August(  )
{
	if( dynamic_instability::decision != 2 )
	{
		cout<<" dynamic_instability::decision != 2 "<<endl;
		cout<<" MatrixXd Cell_two_IS::dynamic_instability_MT_grow_August(  ) "<<endl;
		throw("");
	}



	unsigned int Number_of_MT = sim_of_Cell::Microtubules;
	//If the switch is set to true, the function takes the ammount of microtubules and alternates their lenght
	//with a process simulating the dynamic instability
	//It produces the probability distribution
	if( dynamic_instability::calibration_switch == true  ) 
	{	
		Number_of_MT = 5e5;
	}	
	MatrixXd Matrix_MTs = MatrixXd::Zero( Number_of_MT , 2 );	


	//Initialization of length-set to be the half of the cell circumference
	//All microtubules are growing	
	double initial_lenght = sim_of_Cell::PI * Cell_parametres::A_AXIS;//sim_of_Cell::PI * this->a_axis
	for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
	{
		Matrix_MTs( micro_index , 0 ) = initial_lenght;
		Matrix_MTs( micro_index , 1 ) = 1.0;	
	}

        std::uniform_real_distribution<> distribution{ 0 , 1 };
        unsigned int number_of_generator = omp_get_thread_num();


	//Dynamic instability process
	unsigned int number_of_time_steps = (unsigned int)( dynamic_instability::Time_Of_initial_process / dynamic_instability::Time_Step_Of_in_pr );
	for( unsigned int time_step_id = 0 ; time_step_id <  number_of_time_steps ; time_step_id ++ )
	{

		//grow and shrinkage
		for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
		{
			//Grow		
			if( abs( Matrix_MTs( micro_index , 1 ) - 1.0 ) < 1e-6 )
			{
				Matrix_MTs( micro_index , 0 ) = Matrix_MTs( micro_index , 0 ) + dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::polymerization_constant;
			}
			//Shrinkage			
			else if( abs( Matrix_MTs( micro_index , 1 ) - 2.0 ) < 1e-6 )
		        {
			Matrix_MTs( micro_index , 0 ) = Matrix_MTs( micro_index , 0 ) - dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::shrinking_constant;				
			}
			else
			{
				cout<<"void Cell_two_IS::dynamic_instability_MT_grow(  )"<<endl;
				cout<<"Matrix_MTs( micro_index , 1 ) =  "<<Matrix_MTs( micro_index , 1 )<<endl;
				cout<<"ERROR_ID = 653131313513"<<endl;
				throw("");
			}
			//Too short, growing			
			if( Matrix_MTs( micro_index , 0 )  <=  dynamic_instability::coefficient_4 )
			{
				Matrix_MTs( micro_index , 1 ) = 1.0;
			}


		}		


		//Switching
		for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ;  micro_index ++ )
		{
			//Growing switch to shrinkage
			if( abs( Matrix_MTs( micro_index , 1 ) - 1.0 ) < 1e-6 )
			{
				double Lenght_tmp = Matrix_MTs( micro_index , 0 );
				double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability <  dynamic_instability::rate_based_on_Length_August( Lenght_tmp ) * dynamic_instability::Time_Step_Of_in_pr )
				{
					Matrix_MTs( micro_index , 1 ) = 2;
				}
			}
			//Shrinkage to growing			
			else if( abs( Matrix_MTs( micro_index , 1 ) - 2.0 ) < 1e-6 )
			{
				double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
				if( probability <  dynamic_instability::rescue_rate * dynamic_instability::Time_Step_Of_in_pr )
				{
					Matrix_MTs( micro_index , 1 ) = 1;
				}
			}
		}		
	}

	//Saving the lengths
	if( dynamic_instability::calibration_switch == true  ) 
	{	
        
		char name_of_text_file [100];
		sprintf ( name_of_text_file , "/localdisk/OLD-LD/hornak/matlab_scripts_lenght_distribution/micro_Lenght_distribution_dec_3_2.txt" );
		FILE *out;
		out = fopen( name_of_text_file , "w");

		for( unsigned int row = 0 ; row < Matrix_MTs.rows() ; row ++ )
		{
			fprintf(out,"%20.10f %10.5f\n",Matrix_MTs( row , 0 )  , Matrix_MTs( row ,  1 ) );
		}
		fclose( out );
	}
	//Basic controls
	for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
	{
		double length = Matrix_MTs( micro_index , 0 );
		if( length < dynamic_instability::coefficient_4 - 2.0 * dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::shrinking_constant )
		{
			cout<<" length < dynamic_instability::coefficient_4 - 2.0 * dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::shrinking_constant  "<<endl;
			cout<<"MatrixXd Cell_two_IS::dynamic_instability_MT_grow_August	(  )"<<endl;
			cout<<"ERROR_ID = 946168131614615"<<endl;
			throw("");	
		} 
	}
	//Slight increase
	for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
	{
		double length = Matrix_MTs( micro_index , 0 );
		if( length < dynamic_instability::coefficient_4 )
		{
			cout<<"----------------------------------"<<endl;
			cout<<"dynamic_instability::coefficient_4 = "<<dynamic_instability::coefficient_4<<endl;
			Matrix_MTs( micro_index , 0 ) = dynamic_instability::coefficient_4 + dynamic_instability::polymerization_constant * 0.1; //increase for half a second
		} 
	}
	//The last check, if the microtubule is too short, the error will stop the simulation
	for( unsigned int micro_index = 0 ; micro_index < Number_of_MT ; micro_index ++ )
	{
		double length = Matrix_MTs( micro_index , 0 );
		if( length < dynamic_instability::coefficient_4 )
		{
			cout<<" length < dynamic_instability::coefficient_4 - 2.0 * dynamic_instability::Time_Step_Of_in_pr * dynamic_instability::shrinking_constant  "<<endl;
			cout<<"MatrixXd Cell_two_IS::dynamic_instability_MT_grow_August(  )"<<endl;
			cout<<"ERROR_ID = 1296126849614584941515"<<endl;
			throw("");	
		} 
	}
	return Matrix_MTs;
}




//Implementation of the dynamic instability
void Cell_two_IS::dynamic_instability_August()
{

  double probability_rescue_per_time_step = dynamic_instability::rescue_rate * sim_of_Cell::time_Step;
  double probability_catastrophe_per_time_step = dynamic_instability::catastrophe_rate * sim_of_Cell::time_Step;
  unsigned int number_of_generator = omp_get_thread_num();



  std::uniform_real_distribution<> distribution{ 0 , 1 };
  for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
  {
  	//Only microtubules that are free or attached to cortical sliding dyneins grow - they have a free plus-end  
	if( ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 0 ) && ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 9 ) )
	{
		continue;
	}

	if( ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0 ) || ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 9 ) )
	{
	        //Stable MT will either switch to growing, shrinking 	
		if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 0  )
		{
		    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
			double L_outside_of_MTOC = this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC();
			double probability_catastrophe_depending_on_lenght = dynamic_instability::rate_based_on_Length_August( L_outside_of_MTOC );
			double probability_catastrophe_depending_on_lenght_per_time_step = probability_catastrophe_depending_on_lenght * sim_of_Cell::time_Step;

			if( probability < probability_rescue_per_time_step )
			{
				this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
			}
			else if( ( probability >= probability_rescue_per_time_step ) && ( probability < probability_rescue_per_time_step + probability_catastrophe_depending_on_lenght_per_time_step )  )
			{

				this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
			}
			else
			{

			}
		}
		//Growing microtubule can switch to shrinking		
		else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
		{		
			double L_outside_of_MTOC = this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC();
			double probability_catastrophe_depending_on_lenght = dynamic_instability::rate_based_on_Length_August( L_outside_of_MTOC );
			double probability_catastrophe_depending_on_lenght_per_time_step = probability_catastrophe_depending_on_lenght * sim_of_Cell::time_Step;
		    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
			if( probability < probability_catastrophe_depending_on_lenght_per_time_step )
			{
				this->array_Of_Microtubules[ micro_id ].set_growth_index( 2 );
			}
		}
		//Shrinking MTs can switch to growing		
		else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
		{
		    	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
			if( probability < probability_rescue_per_time_step )
			{
				this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
			}
			if( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() < dynamic_instability::coefficient_4  )
			{
				this->array_Of_Microtubules[ micro_id ].set_growth_index( 1 );
			}
		}
	}
  }

  	//Growing and shrinking
	for( unsigned int micro_id = 0 ; micro_id < this->number_of_microtubules ; micro_id ++ )
	{
		if( ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 0 ) && ( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 9 ) )
		{
			continue;
		}

		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 0  )
		{
			//Growing of microtubules
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{

				this->array_Of_Microtubules[ micro_id ].microtubule_growth_MT_0_9();
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	

			}
			//Shrinking
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{
				this->array_Of_Microtubules[ micro_id ].microtubule_shrinkage_MT_0_9();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	
			}
		}

		//When the microtubule is attached to cortical sliding dyneins, they will detach, if the length is shorter than abscissa
		if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() == 9  )
		{
			if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 1 )
			{
				this->array_Of_Microtubules[ micro_id ].microtubule_growth_MT_0_9();
				this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0();
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	
			}
			else if( this->array_Of_Microtubules[ micro_id ].get_growth_index() == 2 )
			{

				    this->array_Of_Microtubules[ micro_id ].microtubule_shrinkage_MT_0_9();
				    this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				    this->array_Of_Microtubules[ micro_id ].control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9();
				    //This is the point where the dyneins detach				    
				    std::vector< Vector3d > motors_from_micro = this->array_Of_Microtubules[ micro_id ].detachmet_due_to_shrinkage_of_MT_9( );
				    //They have to be returned to either the first, or the second IS
				    //It is similar to the function stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS				    
				    this->array_Of_Microtubules[ micro_id ].set_lenght_of_tangents("void Cell_two_IS::dynamic_instability_linear()");
				    this->array_Of_Microtubules[ micro_id ].set_lenght_of_micro_0_outside_MTOC( this->array_Of_Microtubules[ micro_id ].get_lenght_of_microtubule_outside_MTOC() );	

				    if( motors_from_micro.size() == 0 )
                  		    {
					continue;
				    }	
	
				    std::vector<Vector3d> points_inside;
				    std::vector<Vector3d> returned_outside_IS = this->IS_Cortic_Sl.control_IS_Cortical_Sliding_points_2_b( motors_from_micro , points_inside );
				    std::vector<Vector3d> margin_points;
				    std::vector<Vector3d> points_inside_inner_ring = this->IS_Cortic_Sl.control_of_inside_and_margins( points_inside , margin_points );
				
				    if( points_inside_inner_ring.size() > 0 )
				    {	
				    	this->project_and_add_points_to_surface( points_inside_inner_ring );
				    }	

				    if( margin_points.size() > 0 )
				    {
					//zlabava
					ISCorticalSl IS_tmp_1 = this->get_IS_cortical_sliding_first();
					std::vector<Vector3d> margin_points_2 = this->density_surface_dynein.create_cortical_sliding_points( margin_points.size() , IS_tmp_1 );
					//std::vector<Vector3d> margin_points_2 = this->IS_Cortic_Sl.create_IS_Cortical_Sliding_points( margin_points.size() );
					this->project_and_add_points_to_surface( margin_points_2 );
				    }	

				    std::vector<Vector3d> points_inside_2;
				    std::vector<Vector3d> returned_outside_IS_2 = this->IS_Cortic_Sl_2.control_IS_Cortical_Sliding_points_2_b( returned_outside_IS , points_inside_2 );
				    std::vector<Vector3d> margin_points_IS_2;
				    std::vector<Vector3d> points_inside_inner_ring_2 = this->IS_Cortic_Sl_2.control_of_inside_and_margins( points_inside_2 , margin_points_IS_2 );
				
				    if( points_inside_inner_ring_2.size() > 0 )
				    {	
				    	this->project_and_add_points_to_surface( points_inside_inner_ring_2 );
				    }	

				    if( margin_points_IS_2.size() > 0 )
				    {
					ISCorticalSl IS_tmp_2 = this->get_IS_cortical_sliding_second();      
					std::vector<Vector3d> margin_points_2 = this->density_surface_dynein.create_cortical_sliding_points( margin_points_IS_2.size() , IS_tmp_2 );
					this->project_and_add_points_to_surface( margin_points_2 );
				    }	

				    if( returned_outside_IS_2.size() > 0 )
				    {	
				    	cout<<"void Cell_two_IS::dynamic_instability()"<<endl;
				    	cout<<"returned_outside_IS_2.size() > 0"<<endl;
					throw("");
				    }	
			
			}
		}



	}

}


































//Mid-step algorithm for the cell with two IS	
void Cell_two_IS::MidStep_two_IS()
{

//MICROTUBULES
    //The arrays store the forces, coordinates of all microtubules	
    MatrixXd* random_Forces = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* original_Coordinates = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* external_Forces = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* dynein_Forces = new MatrixXd[ this->number_of_microtubules ];

//MTOC
    //Matrix stores the forces, coordinates of the MTOC
    MatrixXd MTOC_original_coordinates = this->MTOC.get_coordinates();
    MatrixXd MTOC_forces = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_wall_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_nucleus_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );


//Saves original coordinates
for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
    	original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();
}

//Microtubules are unattached, attached to capture-shrikage and cortical sliding
//Since the acting forces are different, algorithm is adapted 
for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
	if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() < 1 )
	{
            //External forces will contain all forces from every structure in cell (dynein, wall, MTOC )
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //Cell wall interactions
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;


            //Nucleus interactions
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC interaction
            //Force acting on the first bead of the microtubule 
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            //Force acting on the second bead of the microtubule             
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );
	    //Calculate forces acting on two first beads of the microtubule and two points(sprouting, rear point)
            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );

            if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2  )
            {
                //I get two points to which the microtubule is connected            
                unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
                unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

                for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
                {
                    //Writing forces of the Microtubule                
                    external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                    if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                    {
                          external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                    }
                }
                //Writing forces of the MTOC                 
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

            //external forces will contain all forces from every structure in cell ( dynein , wall , MTOC )
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //Cell wall 
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

	    //Nucleus
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC
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
	
	    //Dynein forces
            MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
            external_Force = external_Force + force_dynein_IS_surface;

            external_Forces[ microtubule ] = external_Force;

    }

    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
	{
            //external forces will contain all forces from every structure in cell ( dynein , wall , MTOC )
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //Cell wall
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

	    //Nucleus interaction
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC
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



            MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
            dynein_Forces[ microtubule ] = force_dynein_IS_surface;

            external_Force = external_Force + force_dynein_IS_surface;
            external_Forces[ microtubule ] = external_Force;
    }
   else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 )
	{
            //external forces will contain all forces from every structure in cell ( dynein , wall , MTOC )
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //Cell wall
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

	    //Nucleus interaction
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC
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

            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );


            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;


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


            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

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

    }

    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
	{
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

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
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }


            MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );

            external_Forces[ microtubule ] = external_Force;
    }
    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 )
	{
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

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
    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 40 )
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



//Controls the number of dynein in both IS
void Cell_two_IS::test_dynein( string param )
{
	//The number of capture-shrinkage dynein acting on microtubules in the first IS		
	unsigned int number_dynein_20 = this->get_number_of_dynein_motors_capture_shrinkage_1(  );
	//The number of capture-shrinkage dynein in the second IS			
	unsigned int number_dynein_20_IS = this->Capture_Shrinkage_dynein.get_all_dynein_point(  ).size();
	//The original number of dyneins	
	unsigned int basic_number_capt_1 = this->Capture_Shrinkage_dynein.get_original_number_capt();
	
	//The number of capture-shrinkage dynein acting on microtubules in the second IS		
	unsigned int number_dynein_40 = this->get_number_of_dynein_motors_capture_shrinkage_2(  );
	//The number of capture-shrinkage dynein in the second IS			
	unsigned int number_dynein_40_IS = this->Capture_Shrinkage_dynein_2.get_all_dynein_point(  ).size();
	//The original number of dyneins
	unsigned int basic_number_capt_2 = this->Capture_Shrinkage_dynein_2.get_original_number_capt();	
	
	
	//Control of capture-shrinkage dyneins
	if( number_dynein_20 + number_dynein_20_IS != basic_number_capt_1 )
	{
		cout<<param<<endl;
		cout<<"number_dynein_20 = "<<number_dynein_20<<endl;
		cout<<"number_dynein_20_IS = "<<number_dynein_20_IS<<endl;		
		cout<<"basic_number_capt_1 = "<<basic_number_capt_1<<endl;				
		cout<<"number_dynein_20 + number_dynein_20_IS ! = basic_number_capt_1 "<<endl;
		cout<<"number_dynein_40 = "<<number_dynein_40<<endl;
		cout<<"number_dynein_40_IS = "<<number_dynein_40_IS<<endl;		
		cout<<"basic_number_capt_2 = "<<basic_number_capt_2<<endl;			
		cout<<"void Cell_two_IS::test_dynein( string param )"<<endl;
		cout<<"cell_id = "<<this->cell_id<<endl;
		cout<<"ERROR_ID = 364156456132516"<<endl;
		throw("");			
	}

	
	if( number_dynein_40 + number_dynein_40_IS != basic_number_capt_2 )
	{
		cout<<param<<endl;
		cout<<"number_dynein_40 = "<<number_dynein_40<<endl;
		cout<<"number_dynein_40_IS = "<<number_dynein_40_IS<<endl;		
		cout<<"basic_number_capt_2 = "<<basic_number_capt_2<<endl;		
		cout<<" number_dynein_40 + number_dynein_40_IS ! = basic_number_capt_2"<<endl;
		cout<<"number_dynein_20 = "<<number_dynein_20<<endl;
		cout<<"number_dynein_20_IS = "<<number_dynein_20_IS<<endl;		
		cout<<"basic_number_capt_1 = "<<basic_number_capt_1<<endl;				
		cout<<"cell_id = "<<this->cell_id<<endl;		
		cout<<"void Cell_two_IS::test_dynein( string param )"<<endl;
		cout<<"ERROR_ID = 131364561356"<<endl;
		throw("");			
	}
		
	//The same thing done for cortical slidings	
	std::vector<unsigned int> vector_IS = this->density_surface_dynein.get_number_of_point_in_both_IS(  this->IS_Cortic_Sl , this->IS_Cortic_Sl_2 );
			
    	unsigned int IS_cort_1_counter = 0;	
    	unsigned int IS_cort_2_counter = 0;	    
    	for( unsigned int micro_index = 0 ; micro_index < this->number_of_microtubules ; micro_index ++ )
    	{
        	if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 9  )
        	{

	    		std::vector< Vector3d > motors_from_micro = this->array_Of_Microtubules[ micro_index ].get_Dynein_points_without_erasing(); 
	    	
	    		std::vector<Vector3d> points_inside;
            		std::vector<Vector3d> returned_outside_IS = this->IS_Cortic_Sl.control_IS_Cortical_Sliding_points_2_b_for_dynein_testing( motors_from_micro , points_inside );
	    		IS_cort_1_counter = IS_cort_1_counter + points_inside.size(  );	
	

	    		std::vector<Vector3d> points_inside_2;	
            		std::vector<Vector3d> returned_outside_IS_2 = this->IS_Cortic_Sl_2.control_IS_Cortical_Sliding_points_2_b_for_dynein_testing( returned_outside_IS , points_inside_2 );
	    		IS_cort_2_counter = IS_cort_2_counter + points_inside_2.size(  );	
	    
	    		if( points_inside.size(  ) + points_inside_2.size(  ) != motors_from_micro.size() )
	    		{
	    			cout<<"points_inside.size(  ) + points_inside_2.size(  ) != motors_from_micro.size()"<<endl;
	    			cout<<"void Cell_two_IS::test_dynein( string param )"<<endl;
	    			cout<<"ERROR_ID = 6115685154681"<<endl;
	    			throw("");	
	    		}
       		}
     	}	
     	
     	unsigned int original_cort_1 = this->density_surface_dynein.get_original_number_cort_1();
     	unsigned int original_cort_2 = this->density_surface_dynein.get_original_number_cort_2();
     	
	    	
     	if( IS_cort_1_counter + vector_IS[ 0 ] !=  original_cort_1 )
     	{
		cout<<param<<endl;     	
		cout<<"IS_cort_1_counter = "<<IS_cort_1_counter<<endl;
		cout<<"vector_IS[ 0 ] = "<<vector_IS[ 0 ]<<endl;		
		cout<<"original_cort_1 = "<<original_cort_1<<endl;		
		cout<<" IS_cort_1_counter + vector_IS[ 0 ] !=  original_cort_1 "<<endl;
		cout<<"IS_cort_2_counter = "<<IS_cort_1_counter<<endl;
		cout<<"vector_IS[ 1 ] = "<<vector_IS[ 0 ]<<endl;		
		cout<<"original_cort_2 = "<<original_cort_2<<endl;		
		cout<<"cell_id = "<<this->cell_id<<endl;		
		cout<<"void Cell_two_IS::test_dynein( string param )"<<endl;
		cout<<"ERROR_ID = 546413256"<<endl;
		throw("");	     	     	
     	}	
	
	if( IS_cort_2_counter + vector_IS[ 1 ] !=  original_cort_2 )
     	{
		cout<<param<<endl;     	     	
		cout<<"IS_cort_2_counter = "<<IS_cort_1_counter<<endl;
		cout<<"vector_IS[ 1 ] = "<<vector_IS[ 0 ]<<endl;		
		cout<<"original_cort_2 = "<<original_cort_2<<endl;		
		cout<<" IS_cort_2_counter + vector_IS[ 1 ] !=  original_cort_2 "<<endl;
		cout<<"IS_cort_1_counter = "<<IS_cort_1_counter<<endl;
		cout<<"vector_IS[ 0 ] = "<<vector_IS[ 0 ]<<endl;		
		cout<<"original_cort_1 = "<<original_cort_1<<endl;			
		cout<<"void Cell_two_IS::test_dynein( string param )"<<endl;
		cout<<"cell_id = "<<this->cell_id<<endl;		
		cout<<"ERROR_ID = 44645656466"<<endl;
		throw("");	     	     	
     	}	

}







//In this function the cytoskeleton finds a configuration of the smallest energy
//Dynein cannot attach to microtubules
//Microtubules grows and shrinks to assure that their lenghts are distributed according to specified probability density 
//It produces no results
void Cell_two_IS::timeDevelopment_two_Cell_calib( double Time_0 , double Time_1 )
{
        double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++
	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{

       		double time = step * sim_of_Cell::time_Step;
        	Cell_parametres::time = time;
		this->MidStep_two_IS(); //
                this->stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS_calibration_axis_z(  );

		if( dynamic_instability::instability_switch_2 == true )
		{
			if( dynamic_instability::instability_switch_2 == true )
			{
				if( dynamic_instability::decision == 0 )
				{
					this->dynamic_instability();
				}
				else if( dynamic_instability::decision == 1 )
				{
					this->dynamic_instability_linear(  );
				}
			     //The function dynamic_instability_August results in the probability distribution of lenghts described in the publication							
				else if( dynamic_instability::decision == 2 )
				{
					this->dynamic_instability_August(  );
				}
			}
		}
        }
}






void Cell_two_IS::timeDevelopment_Cell_two_IS( double Time_0 , double Time_1 , unsigned int number_of_Pictures )
{
	//Identification of the process
        numerical_result::create_document_file_two_IS();
        
        //Process ID is give. Since it is only one simulation run-0 
	unsigned int key = 0;
        double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++

	unsigned int picture_time_step = 1;

	//Computes when are saved the files for the visualisation of the simulation
	if( numberOfSteps2 > number_of_Pictures )
	{
		picture_time_step = numberOfSteps2 / number_of_Pictures;
	}
	
	
	//Saves the basic shapes of the cell
	this->print_Cell_parametres();
        unsigned int printing_counter = 0;
	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{
		//All objects contain the time information
       		double time = step * sim_of_Cell::time_Step;
		this->MTOC.set_time_clock( time );
		this->set_time_clock( time );

		//Just saving of the files
		////////////////////////////////////////////////////////////////////////////////////////////
		//Creates the files used for the visualisation of the simulation
        	if( numberOfSteps2 >= number_of_Pictures )
		{
			if( step % picture_time_step == 0 )
			{			
				this->print_Cell_two_IS( step / picture_time_step );
			}
		}

		//This shows the progress of the simulation
        	if( step % 1000 == 0)
		{
            		cout<<"key = "<<key<<endl;
            		cout<<"time = "<<time<<endl;
            		//Controls whether NaN is detected. In such a case, the simulation is stopped and the entire configuration of the cell is recorded
            		//Position of every bead and dynein motor is recorded            		
			this->control_NAN( "1000" );
        	}

      		if( step % 200 == 0 ) //50
		{
			//Saves results for statistical analysis		
			this->BIG_PRINT_numerical_two_IS( step , key );			
			//Saves the data about unattached and cortical sliding microtubules			
			this->print_MT_stats( step , key );
        	}
      		if( ( step %  100 ) == 0 ) //14285
		{
			
			//This function assures the uniform distribution of the dyneins in the cortical sliding IS
			this->density_surface_dynein.change_part_of_two_IS_points( this->IS_Cortic_Sl , this->IS_Cortic_Sl_2  );
        	}



      		if( ( step % ( 2 *  7142 ) ) == 0 )
		{
			//Saves the positions of all microtubule beads and all dyneins		
			this->BIG_PRINT_numerical_micro( printing_counter , key );
			//Saves the lenght of all microtubules			
			this->print_MT_stats_lenght( printing_counter , key  );
			printing_counter = printing_counter + 1;
        	}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Microtubules moved under the forces of dynein motors - this->MidStep_two_IS(); 
		//Dyneins are stepping on them - this->stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS();
		//Dyneins attach to microtubules - this->catch_pair_abscissa_real_dynein(); and this->microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro_two_IS();
		//Microtubules grow and shrink - this->dynamic_instability_August(  );
		this->MidStep_two_IS(); //

		//Stepping and detechment of dyneins
                this->stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS();

		//Attachment of cortical sliding dyneins on microtubules
		if( ( sim_of_Cell::density_of_dynein_motor_Cortical_Sliding > 1.0 ) || ( sim_of_Cell::density_of_dynein_motor_Cortical_Sliding_2 > 1.0 )   )
		{
			this->catch_pair_abscissa_real_dynein();
		}

		//Capture-shrinkage
		if( ( sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage > 1.0 ) || ( sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage_2 > 1.0 )  )
		{
			//This controls whether the microtubule attaches to capture-shrinkage dynein		
	       		this->check_if_micro_is_caught_in_IS_two_IS( );
	       		//Attachment of new dyneins on microtubules
			this->microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro_two_IS();
		}
		//Just controls the detachment of capture-shrinkage microtubules		
		this->control_dynein_detachment();
		//Implementation of dynamic instability
		//Three sets of parameters were used		
		if( dynamic_instability::instability_switch_2 == true )
		{
			if( dynamic_instability::instability_switch_2 == true )
			{
				if( dynamic_instability::decision == 0 )
				{
					this->dynamic_instability();
				}
				else if( dynamic_instability::decision == 1 )
				{
					this->dynamic_instability_linear(  );
				}
				//This implementation uses the parameters described in the publication				
				else if( dynamic_instability::decision == 2 )
				{
				        //Implementation of the dynamic instability				
					this->dynamic_instability_August(  );
				}
			}
		}

      		if( ( step %  100 ) == 0 )
		{
			this->test_dynein("0000000000000");
        	}
        }
}








//The same but producing the numerical results
void Cell_two_IS::timeDevelopment_Cell_two_IS_numerical_results( double Time_0 , double Time_1 , unsigned int key )
{
	//Calculates the number of steps
    	double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++
	unsigned int print_integer = 0;
	//This saves the front center of the IS
	this->print_center_numerical_parameters( key );
	unsigned int printing_counter = 0;


	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{
		//Sets the clock
       		double time = step * sim_of_Cell::time_Step;
		this->MTOC.set_time_clock( time );
		this->set_time_clock( time );

		//Controls if NaN is detected
		//If it is, the simulation is stopped
        	if( step % 100000 == 0)
		{
            		cout<<"key = "<<key<<endl;
            		cout<<"time = "<<time<<endl;
			this->control_NAN( std::to_string( key ) );
        	}

      		if( step % 50 == 0 ) //50
		{
			//Saves the numerical results
			this->BIG_PRINT_numerical_two_IS( step , key );
			//this->print_number_of_shorter_microtubules_two_IS( step , key );
        	}


	
      		if( ( step % ( 10000 ) ) == 0 )
		{
			if( sim_of_Cell::writing_switch == true  )
			{	//Saves the entire configuration of the cytoskeleton and dynein motors
				this->BIG_PRINT_numerical_micro( printing_counter , key );
			}	
			//Saves the information about the microtubule lengths, growth and shrinkage
			this->print_MT_stats( printing_counter , key );
			this->print_MT_stats_lenght( printing_counter , key  );
			printing_counter = printing_counter + 1;
        	}


		//Makes sure that dyneins are unifromly distributed 
      		if( ( step %  100 ) == 0 ) //14285
		{
			this->density_surface_dynein.change_part_of_two_IS_points( this->IS_Cortic_Sl , this->IS_Cortic_Sl_2  );
        	}


		//The same as in timeDevelopment_Cell_two_IS, see comments there
		this->MidStep_two_IS(); //

                this->stepping_and_detachment_of_all_microtubule_projection_real_dynein_two_IS();

		if( ( sim_of_Cell::density_of_dynein_motor_Cortical_Sliding > 1.0 ) || ( sim_of_Cell::density_of_dynein_motor_Cortical_Sliding_2 > 1.0 )   )
		{
			this->catch_pair_abscissa_real_dynein();
		}

		if( ( sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage > 1.0 ) || ( sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage_2 > 1.0 )  )
		{
	       		this->check_if_micro_is_caught_in_IS_two_IS( );
	       		
			this->microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro_two_IS();
		}


		//v pohode, zkontrolovany	
		this->control_dynein_detachment();



		if( dynamic_instability::instability_switch_2 == true )
		{
			if( dynamic_instability::instability_switch_2 == true )
			{
				if( dynamic_instability::decision == 0 )
				{
					this->dynamic_instability();
				}
				if( dynamic_instability::decision == 1 )
				{
					//v pohode, zkontrolovany	
					this->dynamic_instability_linear(  );
				}
				if( dynamic_instability::decision == 2 )
				{
					this->dynamic_instability_August(  );
				}

			}
		}

      		if( ( step %  100 ) == 0 )
		{
			this->test_dynein("basic test");
        	}
		


	}

}



































Cell_two_IS::~Cell_two_IS()
{

	delete[] array_Of_Microtubules;
	array_Of_Microtubules = NULL;
	// TODO Auto-generated destructor stub
}
