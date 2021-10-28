#include "numerical_result.h"

namespace numerical_result 
{
    //This is the simulation of one process	
    void repolarization_sample_two_IS( unsigned int index )
    {
    	    //This is the number of microtubules in the cytoskeleton - 100
	    unsigned int number_Of_Microtubules = sim_of_Cell::Microtubules;
	    //This is the number of microtubules sprouting in specified directions from the MTOC
	    //This is for the study of the asymetrical cytoskeleton and so far was not used to generate publicated results 
	    unsigned int number_Of_extra_Microtubules = 0;
	    double Time_0 = 0.0;
	    //The duration of the simulation
	    double Time_1 = sim_of_Cell::Time_1;

	    //Creation of the cell with two IS
	    //The cytosleton is not in the configuration of the lowest energy...
	    //index - Specific id of the process
	    Cell_two_IS cell( number_Of_Microtubules , number_Of_extra_Microtubules , index );
	    //...therefore, 10 seconds are given to find physically realistic confifguration
	    //No data are collected during the timeDevelopment_two_Cell_calib
	    //Microtubules remain free, unatached to dynein
	    cell.timeDevelopment_two_Cell_calib( Time_0 , 10 ); //Time_1
	    //Beginning of the process
	    //Microtubules can attach
	    //Data are collected
	    cell.timeDevelopment_Cell_two_IS_numerical_results( Time_0 , Time_1 , index );
    }	




   //This function only executes all the simulation runs given by repolarization_sample_two_IS
   void test_parallel_two_IS( )
   {
        //Every simulation  run has its specific ID and the results are stored in the files with ID
        //The IDs are chosen randomly in order to avoid the accidental duplication
    	std::uniform_int_distribution<> distribution{ 0 , 10000 };
    	unsigned int number_of_generator = omp_get_thread_num();
	unsigned int lower_boundary = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] ); 
	lower_boundary = lower_boundary * 100;

	//The creation of the file with the parameters of the simulation
	//The data in the file are sufficient to reproduce the results
	numerical_result::create_document_file_two_IS();
	
	//parallel computation of the simulation runs
	#pragma omp parallel for
	for ( int j = lower_boundary ; j < lower_boundary + sim_of_Cell::number_of_samples ; ++j )
	{
		//the variable j is the ID of the simulation run
		repolarization_sample_two_IS( j );
	}
   }


       //Creates the file with the parameters of the simulation 	
       void  create_document_file_two_IS()
	{
		//creates the name of the file
      		std::string name_of_text_file = "./picturesVideos/numerical_results/dokument.txt";
      		ofstream fout;
      		fout.open( name_of_text_file.c_str() , std::fstream::app );
      		//Duration of the simulation      		
      		fout<<sim_of_Cell::Time_1<<endl;  //write to file
      		//time step
      		fout<<sim_of_Cell::time_Step<<endl;  //write to file
      		//Number of microtubules      		
		fout<<sim_of_Cell::Microtubules<<endl;
		//Viscosity
		fout<<sim_of_Cell::viscosity<<endl;
    		fout<<"......................"<<endl;
    		//Density of the cortical sliding dyneins in IS 1
		fout<<sim_of_Cell::density_of_dynein_motor_Cortical_Sliding<<endl;
		//Density of the motors on the plasma membrane
		fout<<sim_of_Cell::density_of_dynein_motor_surface<<endl;
		//Density of the capture-shrinkage dyneinsin IS 1
		fout<<sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage<<endl;
    		fout<<"...."<<endl;
    		//Density of the cortical sliding dyneins in IS 2    		
    		fout<<sim_of_Cell::density_of_dynein_motor_Cortical_Sliding_2<<endl;
		//Density of the capture-shrinkage dyneins in IS 2    		
    		fout<<sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage_2<<endl;
    		fout<<"......................"<<endl;
    		//Radius of the cell
		fout<<Cell_parametres::A_AXIS<<endl;
    		//Radius of the nucleus		
		fout<<Nucleus_parametres::A_AXIS<<endl;
		//Shift of the nucleus on the z-axis: 0 		
		fout<<Nucleus_parametres::z_coordinate<<endl;
                //Number of the original microtubule beads before grow and shrinkage
		fout<<sim_of_Cell::MicrotubulePoints<<endl;
		//Lenght of the segment 
		fout<<sim_of_Cell::resting_distance<<endl;
		//Reduction of the number of beads during the creation
		fout<<sim_of_Cell::REDUCTION<<endl;	
		//Number of the simulation runs
		fout<<sim_of_Cell::number_of_samples<<endl;
		//Radius of the MTOC
                fout<<MTOCparam::MTOC_radius<<endl;
                //Radius of the cortical sliding IS
                fout<<IS_Cortical_Sl_parameter::radius<<endl;	
		fout<<"..........................................."<<endl;
		//Angle between the z-Axis and the axis of the IS1
		fout<<IS_Capture_shrinkage_param::polar_angle<<endl;
		//Angle between the x-Axis and the axis of the IS1		
		fout<<IS_Capture_shrinkage_param::azimutal_angle<<endl;
		//Angle between the z-Axis and the axis of the IS2		
		fout<<IS_Capture_shrinkage_param::polar_angle_2<<endl;
		//Angle between the x-Axis and the axis of the IS2				
		fout<<IS_Capture_shrinkage_param::azimutal_angle_2<<endl;
                fout<<"--------------------------------------------------------------------"<<endl;
		//Parameters of the dynamic instability described in the paper				                
                fout<<"dynamic_Instability"<<endl;	
                fout<<dynamic_instability::coefficient_1<<endl;	
                fout<<dynamic_instability::coefficient_2<<endl;	
                fout<<dynamic_instability::instability_switch<<endl;	
                fout<<dynamic_instability::instability_switch_2<<endl;	
                fout<<dynamic_instability::decision<<endl;	




      		fout.close();
	}




       //The same for the cell with one IS 	
        void  create_document_file()
	{
		//crates the name of the file
      		std::string name_of_text_file = "./picturesVideos/numerical_results/dokument.txt";	
      		ofstream fout;
      		fout.open( name_of_text_file.c_str() , std::fstream::app );	
      		//Duration of the simulation
      		fout<<sim_of_Cell::Time_1<<endl;  //write to file 
      		//Number of microtubules
		fout<<sim_of_Cell::Microtubules<<endl;
		//viscosity
		fout<<sim_of_Cell::viscosity<<endl;
		//Density of the cortical sliding dyneins
		fout<<sim_of_Cell::density_of_dynein_motor_Cortical_Sliding<<endl;
		//Density of dyneins on the surface 
		//It is 0, since the random distribution of dyneins in the surface was not yet studied
		fout<<sim_of_Cell::density_of_dynein_motor_surface<<endl;
		//The density of capture-shrinkage dyneins
		fout<<sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage<<endl;
		//The radius of the plasma membrane
		fout<<Cell_parametres::A_AXIS<<endl;
		//The radius of the nucleus
		fout<<Nucleus_parametres::A_AXIS<<endl;
		//Shift of the nucleus on the z-axis: 0 
		fout<<Nucleus_parametres::z_coordinate<<endl;
		//Angle between z-Axis and the axis of the first IS
		fout<<IS_Capture_shrinkage_param::polar_angle<<endl;
		//Number of the original microtubule beads before grow and shrinkage
		fout<<sim_of_Cell::MicrotubulePoints<<endl;
		//Lenght of the segment 
		fout<<sim_of_Cell::resting_distance<<endl;
		//Reduction of the number of beads during the creation
		fout<<sim_of_Cell::REDUCTION<<endl;
		//Number of the simulation runs
		fout<<sim_of_Cell::number_of_samples<<endl;
		//Radius of the MTOC
                fout<<MTOCparam::MTOC_radius<<endl;
                //Radius of the cortical sliding IS
                fout<<IS_Cortical_Sl_parameter::radius<<endl;	
                fout<<"--------------------------------------------------------------------"<<endl;
                //Parameters of the dynamic instability described in the publication
                fout<<"dynamic_Instability"<<endl;	
                fout<<dynamic_instability::coefficient_1<<endl;	
                fout<<dynamic_instability::coefficient_2<<endl;	
                fout<<dynamic_instability::instability_switch<<endl;	
                fout<<dynamic_instability::instability_switch_2<<endl;	
                //Only the last segment can attach to capture-shrinkage dyneins
                fout<<dynamic_instability::capture_shrinkage_switch<<endl;	
                //the dynamic instability from the publication 
                fout<<dynamic_instability::decision<<endl;	
      		fout.close();
	}	

}



