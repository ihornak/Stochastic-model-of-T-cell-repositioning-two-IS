#include "Surface.h"
#include <stdexcept>

Surface::Surface()
{
//Default constructor, never used
}







//This is the constructor for capture-shrinkage IS
Surface::Surface(  double density , double width , IS_Capt_Shrinkage IS_tmp )
{

    //Just to check IS_tmp if needed
    Vector3d orientation( 0.0 , 0.0 , 0.0 );
    double planar_component = sin( IS_Capture_shrinkage_param::polar_angle );
    double vertical_component = cos( IS_Capture_shrinkage_param::polar_angle );
    double planar_component_x = planar_component * cos( IS_Capture_shrinkage_param::azimutal_angle );
    double planar_component_y = planar_component * sin(IS_Capture_shrinkage_param::azimutal_angle );

    orientation( 0 ) = planar_component_x;
    orientation( 1 ) = planar_component_y;
    orientation( 2 ) = vertical_component;
    orientation = orientation / orientation.norm();

    double minimal_radius = IS_Capture_shrinkage_param::minimal_radius;

    //These are parameters of the map
    //Since the capture-shrinkage has only one key, these parameters do not play a role
    //They has to be defined so that two maps(capture, cortical) have the same structure
    this->A_axis = Cell_parametres::A_AXIS;
    this->B_axis = Cell_parametres::B_AXIS_Lower;
    this->x_width = width;
    this->y_width = width;
    this->z_width = width;
    this->density = density;
    double virtual_sphere_radius = this->A_axis;

    //get Axis from capture shrinkage
    double radius = IS_Capture_shrinkage_param::radius;
    double area = sim_of_Cell::PI * radius * radius;
    double micro_2_to_meter_2 = 1e12;

    //Number of dynein in the IS
    unsigned int number_of_dynein = area * this->density * micro_2_to_meter_2;
    this->original_number_capt = number_of_dynein;
    this->original_number_cort_1 = 0;
    this->original_number_cort_2 = 0;	


    std::vector<Vector3d> vectors_to_be_added = this->create_capture_shrinkage_points( number_of_dynein , IS_tmp );
    this->add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );


}








//Constructs dyneins in two IS
Surface::Surface(  double density , double width , double density_in_IS , ISCorticalSl ISCorticalSl_arg , double density_in_IS_2 , ISCorticalSl ISCorticalSl_arg_2 )
{
    //The parameters of the IS
    unsigned int number_of_generator = omp_get_thread_num();
    this->x_width = width; 
    this->y_width = width;
    this->z_width = width;
    this->density = density;    
    double micro_2_to_meter_2 = 1e12;

     
    //Checking whether plasma membrane does not contain more dyneins than IS
    double difference_density = density_in_IS - density;
    if( difference_density < 0 )
    {
	return;	
    }

    double difference_density_2 = density_in_IS_2 - density;    
    if( difference_density_2 < 0 )
    {
	return;	
    }    
    

    //Creation of points in the first IS
    double area_of_IS_cor_sl = sim_of_Cell::PI * ISCorticalSl_arg.get_radius() * ISCorticalSl_arg.get_radius();
    unsigned int number_of_dynein_IS = ( unsigned int ) ( area_of_IS_cor_sl * difference_density * micro_2_to_meter_2 );
    std::vector<Vector3d>  numbers_in_IS = this->create_cortical_sliding_points(  number_of_dynein_IS  , ISCorticalSl_arg );	
    this->project_and_add_points_to_surface( numbers_in_IS );	
    
    //Creation of points in the second IS
    double area_of_IS_cor_sl_2 = sim_of_Cell::PI * ISCorticalSl_arg_2.get_radius() * ISCorticalSl_arg_2.get_radius();
    unsigned int number_of_dynein_IS_2 = ( unsigned int ) ( area_of_IS_cor_sl_2 * difference_density_2 * micro_2_to_meter_2 );

    //Creation of the cortical sliding dyneins
    std::vector<Vector3d>  numbers_in_IS_2 = this->create_cortical_sliding_points(  number_of_dynein_IS_2  , ISCorticalSl_arg_2 );
    
    //Control, whether the points contain NaN
    for( unsigned int index = 0 ; index < numbers_in_IS_2.size() ; index ++ )
    {	
	/*
	if( std::isnan( numbers_in_IS_2[ index ].norm() )  )
	{
		cout<<" isnan( numbers_in_IS[ index ].norm() ) "<<endl;
		cout<<"Surface::Surface(  double density , double width , double density_in_IS , ISCorticalSl ISCorticalSl_arg , double density_in_IS_2 , ISCorticalSl ISCorticalSl_arg_2 )"<<endl;
		cout<<"ERROR_ID = 615646846564"<<endl;
		throw("");
	}
	*/
    }	
    this->project_and_add_points_to_surface( numbers_in_IS_2 );	
    
    //Saves the number of points for the control
    this->original_number_capt = 0;
    this->original_number_cort_1 = number_of_dynein_IS;
    this->original_number_cort_2 = number_of_dynein_IS_2;	        
}




//Makes sure that the point is located on plasma membrane
Vector3d Surface::project_point_on_surface( Vector3d position )
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


////Adds points to the surface
void Surface::project_and_add_points_to_surface( std::vector< Vector3d > points )
{

    Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );
    for( unsigned int point_id = 0 ; point_id < points.size() ; point_id ++ )
    {
	Vector3d position_updated = points[ point_id ];
        Vector3d divide_position_updated = position_updated + add_to_get_segment;
        //cout<<position_updated<<endl;

        unsigned int x_index = divide_position_updated( 0 ) / this->x_width;
        unsigned int y_index = divide_position_updated( 1 ) / this->y_width;
        unsigned int z_index = divide_position_updated( 2 ) / this->z_width;

        unsigned int ID_map_index = x_index * this->neccessary_dimension * this->neccessary_dimension + y_index * this->neccessary_dimension + z_index;
        this->surface[ ID_map_index ].push_back( position_updated );

    }




}




//Constructor for cortical sliding dyneins on the surface of the cell and in the IS 
Surface::Surface(  double density , double width , string cortical_Sliding , double density_in_IS , ISCorticalSl ISCorticalSl_arg )
{
    //density - density of dynein outside IS
    //density_in_IS - density_of dynein in the IS
    // So far, only density = 0 was used to get statistical results, no dyneins are placed outside the IS 
    //Parameters of the IS

/*
    if( sim_of_Cell::density_of_dynein_motor_Cortical_Sliding < sim_of_Cell::density_of_dynein_motor_surface )
    {
	cout<<"sim_of_Cell::density_of_dynein_motor_Cortical_Sliding < sim_of_Cell::density_of_dynein_motor_surface"<<endl;
        cout<<"ERROR_ID = 65415615616"<<endl;
        cout<<"Surface::Surface(  double density , double width , string cortical_Sliding , double density_in_IS , ISCorticalSl ISCorticalSl_arg )"<<endl;
        throw("");
    }
*/
    this->A_axis = Cell_parametres::A_AXIS;
    this->B_axis = Cell_parametres::B_AXIS_Lower;
    this->x_width = width; 
    this->y_width = width;
    this->z_width = width;
    this->density = density;
    double virtual_sphere_radius = ( this->A_axis + this->B_axis ) / 2.0;
    double surface = 4.0 * sim_of_Cell::PI * virtual_sphere_radius * virtual_sphere_radius;

    double micro_2_to_meter_2 = 1e12;
    unsigned int number_of_dynein = surface * this->density * micro_2_to_meter_2;
    unsigned int counter = 0;

    //Creating points on the surface of the sphere
    std::uniform_real_distribution<> distribution{  ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius };
    unsigned int number_of_generator = omp_get_thread_num();

    std::vector<Vector3d>  numbers_of_vectors_on_surface;
    while( counter < number_of_dynein )
    {

        double x = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        double y = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        double z = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        Vector3d position( x , y , z );

        if( position.norm() > virtual_sphere_radius )
        {
            continue;
        }

        position = position / position.norm() * virtual_sphere_radius;
	numbers_of_vectors_on_surface.push_back( position );
	counter++;

    }
    
    //Add them on the surface of the sphere    
    this->project_and_add_points_to_surface( numbers_of_vectors_on_surface );

    //Filling of IS_cortical sliding
    
    double difference_density = density_in_IS - density;
    Vector3d orientation_of_IS_cortical_sliding = IS_Cortical_Sl_parameter::orientation;
    if( difference_density <= 0 )
    {
	return;
    }
    double area_of_IS_cor_sl = sim_of_Cell::PI * IS_Cortical_Sl_parameter::radius * IS_Cortical_Sl_parameter::radius;
    unsigned int number_of_dynein_IS = ( unsigned int ) ( area_of_IS_cor_sl * difference_density * micro_2_to_meter_2 );

    //Dyneins are created in the IS
    std::vector<Vector3d>  numbers_in_IS = this->create_cortical_sliding_points(  number_of_dynein_IS  , ISCorticalSl_arg );
    //and projected on the plasma membrane    
    this->project_and_add_points_to_surface( numbers_in_IS );
}









Surface::Surface(const Surface& other)
{
    this->A_axis = other.A_axis;
    this->B_axis = other.B_axis;
    this->x_width = other.x_width;
    this->y_width = other.y_width;
    this->z_width = other.z_width;
    this->surface = other.surface;
    this->original_number_capt = other.original_number_capt;
    this->original_number_cort_1 = other.original_number_cort_1;    
    this->original_number_cort_2 = other.original_number_cort_2;  
}

Surface& Surface::operator=( const Surface &tmp )
{
    this->A_axis = tmp.A_axis;
    this->B_axis = tmp.B_axis;
    this->x_width = tmp.x_width;
    this->y_width = tmp.y_width;
    this->z_width = tmp.z_width;
    this->surface = tmp.surface;
    this->original_number_capt = tmp.original_number_capt;
    this->original_number_cort_1 = tmp.original_number_cort_1;    
    this->original_number_cort_2 = tmp.original_number_cort_2;        
    return *this;
}

//Returns the key of the compartment containing the point projected on the plasma membrane
unsigned int Surface::get_dynein_compartment_id_projected( Vector3d point_position )
{
        double BB;
        double B_multiplicator;
        if( point_position( 2 ) < 0 )
        {
            BB = Cell_parametres::B_AXIS_Lower;
            B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        }
        else
        {
            BB = Cell_parametres::B_AXIS;
            B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        }
        double numerator = Cell_parametres::A_AXIS * BB;
        double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;

        double denominator = ( point_position( 0 ) * point_position( 0 ) + point_position( 1 ) * point_position( 1 ) ) * B_multiplicator;
        denominator = denominator + point_position( 2 ) * point_position( 2 ) * A_multiplicator;
        denominator = sqrt( denominator );
        double c_multiplicator = numerator / denominator;

        Vector3d position_updated = point_position * c_multiplicator;
        Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );

        Vector3d position_to_get_index = position_updated + add_to_get_segment;

        unsigned int x_index = position_to_get_index( 0 ) / this->x_width;
        unsigned int y_index = position_to_get_index( 1 ) / this->y_width;
        unsigned int z_index = position_to_get_index( 2 ) / this->z_width;


        unsigned int ID_map_index = x_index * this->neccessary_dimension * this->neccessary_dimension + y_index * this->neccessary_dimension + z_index;
        if( this->surface.count( ID_map_index ) == 0 )
        {
            ID_map_index = 0;
        }

        return ID_map_index;

}







//////////////////////////////////////////////////////////////////////////////////////////////
//It controls, whether the point is contained in the map
//If so, it is erased
bool Surface::dynein_surface_catch_control( Vector3d point_position , Vector3d& Catching_point )
{
    //First, I have to get the key of the compartment that can ontain the point
    if( point_position( 2 ) < -4.5e-6 )
    {
        cout<<"zlabava"<<endl;
        return false;
    }



        double BB;
        double B_multiplicator;
        if( point_position( 2 ) < 0 )
        {
            BB = Cell_parametres::B_AXIS_Lower;
            B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        }
        else
        {
            BB = Cell_parametres::B_AXIS;
            B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        }
        double numerator = Cell_parametres::A_AXIS * BB;
        double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;

        double denominator = ( point_position( 0 ) * point_position( 0 ) + point_position( 1 ) * point_position( 1 ) ) * B_multiplicator;
        denominator = denominator + point_position( 2 ) * point_position( 2 ) * A_multiplicator;
        denominator = sqrt( denominator );
        double c_multiplicator = numerator / denominator;

        Vector3d position_updated = point_position * c_multiplicator;
        //I want to get the right segment
        Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );

        Vector3d position_to_get_index = position_updated + add_to_get_segment;
        unsigned int x_index = position_to_get_index( 0 ) / this->x_width;
        unsigned int y_index = position_to_get_index( 1 ) / this->y_width;
        unsigned int z_index = position_to_get_index( 2 ) / this->z_width;
	//It gets the index of the compartment
        unsigned int ID_map_index = x_index * this->neccessary_dimension * this->neccessary_dimension + y_index * this->neccessary_dimension + z_index;

        try
        {
            std::vector<Vector3d> dynein_on_surface_points = this->surface.at( ID_map_index );
            for( unsigned int point_id = 0 ; point_id < dynein_on_surface_points.size() ; point_id ++ )
            {
                Vector3d dynein_position = dynein_on_surface_points.at( point_id );
		//The points are interchangeable 
                if( ( dynein_position - position_updated ).norm() < IS_Dynein_Cell_surface::dynein_on_surface_catch_radius )
                {
                    Catching_point = dynein_position;
                    //Point is erased  
                    this->erase_dynein_point( ID_map_index , point_id );
                    return true;
                }
            }

        }
        catch (const std::out_of_range oor)
        {
            Catching_point = Vector3d( 0.0 , 0.0 , 0.0 );
            return false;
        }
        return false;
}













//Erase all points in the compartment with the key
void Surface::erase_dynein_point( unsigned int map_key , unsigned int array_index )
{
    //If the map does not contain the compartment, error ends simulation
    try
    {
        this->surface.at( map_key );
    }
    catch (const std::out_of_range& oor)
    {
        std::cerr << "this->surface.at( map_key ) "<< '\n';
        cout<<" wrong key in the c++ map "<<endl;
        cout<<"Surface ERROR_ID = 4843486464838"<<endl;
        throw("");
    }

    if( array_index >= this->surface.at( map_key ).size() )
    {
        cout<<"array_index <= this->surface.at( map_key ).size() in erase_dynein_point( unsigned int map_key , unsigned int array_index )"<<endl;
        cout<<"map_key = "<<map_key<<endl;
        cout<<"this->surface.at( map_key ).size() = "<<this->surface.at( map_key ).size()<<endl;
        cout<<"array_index = "<<array_index<<endl;
        cout<<"Surface ERROR_ID = 68441489348943"<<endl;
        throw("");
    }

    this->surface.at( map_key ).erase( this->surface.at( map_key ).begin() + array_index );

}

void Surface::erase_vector_points_with_key( unsigned int map_key )
{
    try
    {
        this->surface.at( map_key );
    }
    catch (const std::out_of_range& oor)
    {
        std::cerr << "this->surface.at( map_key ) "<< '\n';
        cout<<" wrong key in the c++ map "<<endl;
        cout<<"erase_vector_points_with_key()"<<endl;
        cout<<"Surface ERROR_ID = 46464838"<<endl;
        throw("");
    }
     this->surface[ map_key ].clear();

}


//Get all dyneins in the compartment with the key
//If does not contain any, returns empty vector
std::vector<Vector3d> Surface::get_dynein_points( unsigned int key )
{

    try
    {
        return this->surface.at( key );
    }
    catch( const std::out_of_range& oor )
    {
        //cout<<"std::vector<Vector3d> Surface::get_dynein_points( unsigned int key )"<<endl;
        //cout<<"Surface ERROR_ID = 684735168543483754"<<endl;
        //throw("");
    }
    std::vector<Vector3d> tmp;
    return tmp;

}




//Inserts the vector of points in the compartment with the key
void Surface::set_dynein_points( unsigned int key , std::vector<Vector3d> vectors_arg )
{
    if( vectors_arg.size() == 0 )
    {
	this->surface[ key ].clear();
        return;
    }
    if( vectors_arg.size() > 0 )
    {
        bool answer = this->surface.count( key );
        if( answer == 0 )
        {
            cout<<"You are trying to insert std::vector<Vector3d> insto the map and key does not exist"<<endl;
            cout<<"You are inventing the news segment"<<endl;
            cout<<"map_key = "<<key<<endl;
            cout<<"vectors_arg.size() = "<<vectors_arg.size()<<endl;
            cout<<"Surface ERROR_ID = 6847396434654875"<<endl;
            throw("");
        }
        this->surface[ key ].clear();;
        this->surface[ key ] = vectors_arg;
    }

}

void Surface::add_dynein_point( unsigned int key , Vector3d vector_arg )
{
    /*
    bool answer = this->surface.count( key );
    if( answer == 0 )
    {
        cout<<"Surface::add_dynein_point( unsigned int key , Vector3d vector_arg )"<<endl;
        cout<<"You are trying to insert std::vector<Vector3d> into the map and key does not exist"<<endl;
        cout<<"You are inventing the news segment"<<endl;
        cout<<"map_key = "<<key<<endl;
        unsigned int ERROR_ID = 96816948;
        cout<<"Surface ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }
    */
    this->surface[ key ].push_back( vector_arg );

}

//Gets the number of dyneins in the map
unsigned int Surface::get_dynein_motors_number()
{
    unsigned int sum = 0;
    for (std::map< unsigned int , std::vector<Vector3d> >::iterator it = this->surface.begin(); it!=this->surface.end(); ++it)
    {
        //Iterates through the compartments and adds the number of points    
        std::vector<Vector3d> Vector_of_Points = it->second;
        sum = sum + Vector_of_Points.size();
    }
    return sum;
}

//Creates the dyneins for the capture-shrinkage IS
std::vector<Vector3d> Surface::create_capture_shrinkage_points(  unsigned int number_of_points  , IS_Capt_Shrinkage IS_tmp )
{
    //Points are created close to the cell membrane around the z-axis
    //Then, the points are rotated 
    double minimal_radius = IS_Capture_shrinkage_param::minimal_radius;


    this->A_axis = Cell_parametres::A_AXIS;
    this->B_axis = Cell_parametres::B_AXIS_Lower;

    unsigned int number_of_dynein = number_of_points;
    double virtual_sphere_radius = this->A_axis;

    std::vector<Vector3d> returned_vectors;
    Vector3d z_axis( 0.0 , 0.0 , 1.0 );
    Vector3d axis_IS = IS_tmp.get_axis_of_IS() * ( -1.0 );

    //Computes the angle between the axis of the IS and z-Axis
    double arcosinus = z_axis.dot( axis_IS ) / ( axis_IS.norm() * z_axis.norm() );
	if( abs( arcosinus ) > 1 )
	{
		if( arcosinus < 0 )
		{
			arcosinus = -1.0;
		}
		else
		{
			arcosinus = 1.0;
		}
	}
        double rotation_angle = acos( arcosinus  );
        //I get the axis of rotation        
	Vector3d axis_of_rotation( 0.0 , 0.0 , 0.0  );


	if( abs( rotation_angle - sim_of_Cell::PI ) < 0.1 )
	{
		axis_of_rotation( 0 ) = 1.0;

	}
	else if( abs( rotation_angle ) < 0.01 )
	{
		axis_of_rotation( 0 ) = 1.0;
	}
	else
	{
		axis_of_rotation = z_axis.cross( axis_IS );
        	axis_of_rotation = axis_of_rotation / axis_of_rotation.norm();
	}


   std::uniform_real_distribution<> distribution( 0.0 , 2.0 * sim_of_Cell::PI );
   std::uniform_real_distribution<> distribution_2( minimal_radius , Cell_parametres::A_AXIS );
   unsigned int number_of_generator = omp_get_thread_num();
   //Quaternion for the rotation
    Quaternion<double> q;
    q = AngleAxis<double>( rotation_angle , axis_of_rotation );

    for( unsigned int counter = 0 ; counter < number_of_points ; counter ++  )
    {

    	//The position of the point is determined by the radius_tmp - the distance from the axis of the IS
        double azimutal_angle_tmp = distribution(  mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        double radius_tmp = this->triangular_distribution(  IS_Capture_shrinkage_param::radius );
	double x = radius_tmp * cos( azimutal_angle_tmp );
        double y = radius_tmp * sin( azimutal_angle_tmp );
        Vector3d orientation( 0.0 , 0.0 , 0 );
	orientation( 0 ) = x;
	orientation( 1 ) = y;
        double z = sqrt(  Cell_parametres::B_AXIS *  Cell_parametres::B_AXIS - ( x * x + y * y ) );
        orientation( 2 ) = z;
      	orientation = orientation / orientation.norm();
	      //double random_number = rand_x( minimal_radius , Cell_parametres::A_AXIS );
	double random_number = distribution_2(  mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ]);
	Vector3d point = orientation * random_number;

	//The point is rotated to the axis of the IS
	Vector3d position_updated = q * point;
	returned_vectors.push_back( position_updated );
    }
    return returned_vectors;

}

//Creates cortical sliding points
std::vector<Vector3d> Surface::create_cortical_sliding_points(  unsigned int number_of_points  , ISCorticalSl IS_tmp )
{
	//The cortical sliding dyneins are generated around the z-Axis and than rotated towards the IS
	//The angle is computed
	Vector3d z_axis( 0.0 , 0.0 , 1 );
        Vector3d axis_IS = IS_tmp.get_axis();
	double arcosinus = z_axis.dot( axis_IS ) / ( axis_IS.norm() * z_axis.norm() );
	if( abs( arcosinus ) > 1 )
	{
		if( arcosinus < 0 )
		{
			arcosinus = -1.0;
		}
		else
		{
			arcosinus = 1.0;
		}
	}
        double angle = acos( arcosinus );

	//Axis of rotation is calculated
	Vector3d axis_of_rotation( 0.0 , 0.0 , 0.0  );
	if( abs( angle - sim_of_Cell::PI ) < 0.1 )
	{
		axis_of_rotation( 0 ) = 1.0;

	}
	else if( abs( angle ) < 0.01 )
	{
		axis_of_rotation( 0 ) = 1.0;
	}
	else
	{
		axis_of_rotation = z_axis.cross( axis_IS );
        	axis_of_rotation = axis_of_rotation / axis_of_rotation.norm();
	}

	std::vector<Vector3d> returned_points;
	//Quaternion rotates the points
        Quaternion<double> q;
	q = AngleAxis<double>( angle , axis_of_rotation );
	//cout<<q<<endl;
	std::uniform_real_distribution<> distribution{ 0.0 , 2.0 * sim_of_Cell::PI };
        unsigned int number_of_generator = omp_get_thread_num();
	for( unsigned int number = 0 ; number < number_of_points ; number ++ )
	{
		//The point is determined by the radius from the axis and the angle	
		double azimuthal_angle = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
    		double radius_tmp = this->triangular_distribution(  IS_Cortical_Sl_parameter::radius );
		double x = radius_tmp * cos( azimuthal_angle );
    		double y = radius_tmp * sin( azimuthal_angle );
    		double z = sqrt(  Cell_parametres::B_AXIS *  Cell_parametres::B_AXIS - ( x * x + y * y ) );

		Vector3d bod( x , y , z );
		bod = bod / bod.norm() * Cell_parametres::B_AXIS;
		//Point is rotated		
		Vector3d rotated_point = q * bod;
		returned_points.push_back( rotated_point );
	}
	return returned_points;

}






//Gets all dynein motors from the map
std::vector<Vector3d> Surface::get_all_dynein_point()
{
    std::vector<Vector3d> all_points;
    //For loop iterates throuh the map and dynein motors are added to std::vector<Vector3d> all_points
    for (std::map< unsigned int , std::vector<Vector3d> >::iterator it = this->surface.begin(); it!=this->surface.end(); ++it)
    {
        std::vector<Vector3d> Vector_of_Points = it->second;

        for( unsigned int point_i = 0 ; point_i < Vector_of_Points.size() ; point_i ++ )
        {
    		all_points.push_back( Vector_of_Points[ point_i ] );
        }
    }
    return all_points;
}






//All the points from the map are erased
void Surface::erase_map( )
{
    for (std::map< unsigned int , std::vector<Vector3d> >::iterator it = this->surface.begin(); it!=this->surface.end(); ++it)
    {
	this->erase_vector_points_with_key( it->first );
    }
}



void Surface::change_part_of_IS_points(  ISCorticalSl IS_tmp )
{

        //Procentage of points that will be replaced
	double procentage = sim_of_Cell::surface_replacement_procentage;
	//Gets the points of the map	
	std::vector<Vector3d> dynein_points = this->get_all_dynein_point();
	unsigned int number_of_points = dynein_points.size();
	//The number of replaced motors	
	unsigned int number_of_points_to_be_replaced = procentage * number_of_points;
	if ( number_of_points_to_be_replaced == 0 )
	{
		return;
	}
	//All points in the map are erased 	
	this->erase_map();



	std::vector<Vector3d> new_points = create_cortical_sliding_points( number_of_points_to_be_replaced , IS_tmp );


	int meze = number_of_points - 1;
	std::uniform_int_distribution<> distribution{ 0 ,  meze };
	unsigned int number_of_generator = omp_get_thread_num();
	//Random points are replaced
	for( unsigned int index = 0 ; index < number_of_points_to_be_replaced ; index ++ )
	{
		unsigned int random_int = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
		dynein_points[ random_int ] = new_points[ index ];
	}
	//And added to the surface
	this->project_and_add_points_to_surface( dynein_points );
}

















void Surface::change_part_of_points_surface_and_IS( ISCorticalSl IS_tmp )
{
	unsigned int number_of_generator = omp_get_thread_num();
	double procentage = sim_of_Cell::surface_replacement_procentage;

	std::vector<Vector3d> dynein_points = this->get_all_dynein_point();
	//unsigned int number_of_points = dynein_points.size();

	std::vector<Vector3d> points_outside;
	std::vector<Vector3d> points_inside;
	for( unsigned int index = 0 ; index < dynein_points.size() ; index ++ )
	{
		Vector3d dynein_position = dynein_points[ index ];
		bool answer = IS_tmp.control_IS_Cortical_Sliding_points( dynein_position );
		if ( answer == true )
		{
			points_inside.push_back( dynein_position );
		}
		else
		{
			points_outside.push_back( dynein_position );
		}

	}
	unsigned int number_of_points_in_IS_to_replace = procentage * points_inside.size();
	unsigned int number_of_points_outside_IS_to_replace = procentage * points_outside.size();

	std::vector<Vector3d> new_points_in_IS = create_cortical_sliding_points( number_of_points_in_IS_to_replace , IS_tmp );
	std::vector<Vector3d> new_points_outside_IS = create_Cortical_Sliding_points_outside_IS( number_of_points_outside_IS_to_replace , IS_tmp );

	int meze = points_inside.size() - 1;
	std::uniform_int_distribution<> distribution( 0 , meze );

	for( unsigned int index = 0 ; index <  new_points_in_IS.size() ; index ++ )
	{
		unsigned int random_int = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
		points_inside[ random_int ] = new_points_in_IS[ index ];
	}

	int meze_2 = points_outside.size() - 1;
	std::uniform_int_distribution<> distribution_2( 0 , meze_2 );

	for( unsigned int index = 0 ; index <  new_points_outside_IS.size() ; index ++ )
	{
		unsigned int random_int = distribution_2( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
		points_outside[ random_int ] = new_points_outside_IS[ index ];
	}

	this->erase_map();

	this->project_and_add_points_to_surface( points_inside );
	this->project_and_add_points_to_surface( points_outside );

}











//It does what it says
void Surface::add_dynein_points_to_compartment_with_number( std::vector<Vector3d> vectors  , unsigned int compartment_id )
{
    for( unsigned int counter = 0 ; counter < vectors.size() ; counter ++ )
    {
        Vector3d tmp = vectors[ counter ];
        try
        {
            this->add_dynein_point( compartment_id , tmp );
        }
        catch( unsigned int error_id )
        {
            cout<<"exception"<<endl;
            cout<<"void Surface::add_dynein_points_to_compartment_with_number( std::vector<Vector3d> vectors  , unsigned int compartment_id )"<<endl;
        }
    }

}




//It prints the map to the terminal-just for very basic checking
void Surface::print_inner_map()
{
    for (std::map< unsigned int , std::vector<Vector3d> >::iterator it = this->surface.begin(); it!=this->surface.end(); ++it)
    {
        unsigned int area_key = it->first;
        std::vector<Vector3d> Vector_of_Points = it->second;
        cout<<"................................................."<<endl;
        cout<<" key of area = "<<area_key<<endl;

        for( unsigned int point_i = 0 ; point_i < Vector_of_Points.size() ; point_i ++ )
        {
            cout<<"......."<<endl;
            cout<<Vector_of_Points.at( point_i )<<endl;
        }
        cout<<endl;

    }
}



std::vector<Vector3d> Surface::create_Cortical_Sliding_points_outside_IS( unsigned int number_of_points , ISCorticalSl ISCorticalSl_arg )
{
    unsigned int counter = 0;
    std::vector<Vector3d> returned_vectors;
    double virtual_sphere_radius = ( this->A_axis + this->B_axis ) / 2.0;

    std::uniform_real_distribution<> distribution{  ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius };
    unsigned int number_of_generator = omp_get_thread_num();
    //double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

    while( counter < number_of_points )
    {
        //double x = randx( ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius );
        //double y = randx( ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius );
        //double z = randx( ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius );

        double x = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        double y = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        double z = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        Vector3d position( x , y , z );

        if( position.norm() > virtual_sphere_radius )
        {
            continue;
        }

        double BB;
        double B_multiplicator;
        if( z  < 0 )
        {
            BB = Cell_parametres::B_AXIS_Lower;
            B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        }
        else
        {
            BB = Cell_parametres::B_AXIS;
            B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        }
        double numerator = Cell_parametres::A_AXIS * BB;
        double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;

        double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;
        denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
        denominator = sqrt( denominator );
        double c_multiplicator = numerator / denominator;
        Vector3d position_updated = position * c_multiplicator;
	if( ISCorticalSl_arg.control_IS_Cortical_Sliding_points( position_updated ) == true )
	{
		continue;
	}
	else
	{
                returned_vectors.push_back( position_updated );
		counter++;
	}


    }

    return returned_vectors;

}







void Surface::change_part_of_two_IS_points(  ISCorticalSl IS_tmp_1 , ISCorticalSl IS_tmp_2 )
{
	unsigned int number_of_generator = omp_get_thread_num();
	//The procentage of points, that  will be replaced	
	double procentage = sim_of_Cell::surface_replacement_procentage;
	//Returns all dynein motors from the map	
	std::vector<Vector3d> dynein_points = this->get_all_dynein_point();

	std::vector<Vector3d> points_first_IS;
	std::vector<Vector3d> points_second_IS;
	std::vector<Vector3d> undetermined;
	
	//Decision, whether points belong to the first, or the second IS	
	for( unsigned int index = 0 ; index < dynein_points.size() ; index ++ )
	{
		Vector3d dynein_position = dynein_points[ index ];
		bool answer_1 = IS_tmp_1.control_IS_Cortical_Sliding_points( dynein_position );
		bool answer_2 = IS_tmp_2.control_IS_Cortical_Sliding_points( dynein_position );
		if ( answer_1 == true )
		{
			points_first_IS.push_back( dynein_position );
		}
		else if ( answer_2 == true )
		{
			points_second_IS.push_back( dynein_position );

		}
		else
		{
			//Since the smallest angle between two IS in our simulation is pi/2, the IS are always sufficiently appart to distinguish the points
			//In such a case, undetermined is unused		
			undetermined.push_back( dynein_position );
		}

	}

	//The number of dynein that will be replaced in the first and the second IS
	unsigned int number_of_IS_1111_to_replace = procentage * points_first_IS.size();
	unsigned int number_of_IS_2222_to_replace = procentage * points_second_IS.size();

	//New points are created 
	std::vector<Vector3d> new_points_in_IS_1111 = this->create_cortical_sliding_points( number_of_IS_1111_to_replace , IS_tmp_1 );
	std::vector<Vector3d> new_points_in_IS_2222 = this->create_cortical_sliding_points( number_of_IS_2222_to_replace , IS_tmp_2 );

        //The points in the first IS are created
	if( points_first_IS.size() > 0 )
	{
		int meze_1 = points_first_IS.size() - 1;
		std::uniform_int_distribution<> distribution_1( 0 , meze_1 );

		//menim body v prvni IS
		for( unsigned int index = 0 ; index <  new_points_in_IS_1111.size() ; index ++ )
		{
			unsigned int random_int = distribution_1(mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ]);
			points_first_IS[ random_int ] = new_points_in_IS_1111[ index ];
		}
	}
	//The points in the second IS are created
	if( points_second_IS.size() > 0 )
	{
		int meze_2 = points_second_IS.size() - 1;
		std::uniform_int_distribution<> distribution_2( 0 , meze_2 );

		//Random points in the second IS are changed
		for( unsigned int index = 0 ; index <  new_points_in_IS_2222.size() ; index ++ )
		{
			unsigned int random_int = distribution_2( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
			points_second_IS[ random_int ] = new_points_in_IS_2222[ index ];
		}
	}

	//If necessary, it controls nans
	if( points_first_IS.size() >=1  )
	{
		for( unsigned int index = 0 ; index < points_first_IS.size() ; index ++ )
		{
			Vector3d point_tmp = points_first_IS[ index ];
/*
			if( isnan(point_tmp.norm() ) )
			{
				cout<<"------------------------------------------------------"<<endl;
			}
*/
		}
	}
	//The points in the second IS are created
	if( points_second_IS.size() >=1  )
	{
		for( unsigned int index = 0 ; index < points_second_IS.size() ; index ++ )
		{
/*
			Vector3d point_tmp = points_second_IS[ index ];
			if( isnan(point_tmp.norm() ) )
			{
				cout<<"------------------------------------------------------"<<endl;
			}
*/
		}
	}

	//The original map is erased
	this->erase_map();
	//Points are projected and added
	this->project_and_add_points_to_surface( points_first_IS );
	this->project_and_add_points_to_surface( points_second_IS );
}





//Gets the number of dynein points in both IS - just for control
std::vector<unsigned int> Surface::get_number_of_point_in_both_IS(  ISCorticalSl IS_tmp_1 , ISCorticalSl IS_tmp_2 )
{
	std::vector<Vector3d> points_first_IS;
	std::vector<Vector3d> points_second_IS;
	std::vector<Vector3d> undetermined;
	
	//Returns all the points in the map
	std::vector<Vector3d> dynein_points = this->get_all_dynein_point();
	for( unsigned int index = 0 ; index < dynein_points.size() ; index ++ )
	{
		Vector3d dynein_position = dynein_points[ index ];
		//Controls to which IS dynein belongs to 
		bool answer_1 = IS_tmp_1.control_IS_Cortical_Sliding_points( dynein_position );
		bool answer_2 = IS_tmp_2.control_IS_Cortical_Sliding_points( dynein_position );
		if ( answer_1 == true )
		{
			points_first_IS.push_back( dynein_position );
		}
		else if ( answer_2 == true )
		{
			points_second_IS.push_back( dynein_position );

		}
		else
		{
			undetermined.push_back( dynein_position );
		}
	}
	
	//If some point is located outside IS, the simulation is interupted
	//This is used only when no point is generated on the membrane
	if( undetermined.size() != 0 )
	{
		cout<<"undetermined.push_back( dynein_position );"<<endl;
		cout<<"void Surface::get_number_of_point_in_both_IS(  ISCorticalSl IS_tmp_1 , ISCorticalSl IS_tmp_2 )"<<endl;
		cout<<"ERROR_ID = 1633131511313"<<endl;
		throw("");		
	}
	
	std::vector<unsigned int> results;
	results.push_back( points_first_IS.size() );
	results.push_back( points_second_IS.size() );	
	return results; 
}










//Generate the number from the triangular distribution between zero and upper boundary
double Surface::triangular_distribution( double upper_boundary )
{

  std::uniform_real_distribution<> distribution{ 0 , 1 };
  unsigned int number_of_generator = omp_get_thread_num();
  double R = upper_boundary;
  double x = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

  double result = sqrt( x ) * R;
  return result;
}





Surface::~Surface()
{

}




//Constructor for cortical sliding IS just for the purposes of testing
Surface::Surface(  double density , double width , string cortical_Sliding , double density_in_IS )
{

    //Basic parameters of the cortical sliding IS 
    this->A_axis = Cell_parametres::A_AXIS;
    this->B_axis = Cell_parametres::B_AXIS_Lower;
    this->x_width = width; 
    this->y_width = width;
    this->z_width = width;
    this->density = density;
    double virtual_sphere_radius = ( this->A_axis + this->B_axis ) / 2.0;
    double surface = 4.0 * sim_of_Cell::PI * virtual_sphere_radius * virtual_sphere_radius;

    double micro_2_to_meter_2 = 1e12;
    unsigned int number_of_dynein = surface * this->density * micro_2_to_meter_2;
    unsigned int counter = 0;


    std::uniform_real_distribution<> distribution{  ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius };
    unsigned int number_of_generator = omp_get_thread_num();
    while( counter < number_of_dynein )
    {
        double x = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        double y = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        double z = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
	Vector3d position( x , y , z );

        if( position.norm() > virtual_sphere_radius )
        {
            continue;
        }


        double BB;
        double B_multiplicator;
        if( z  < 0 )
        {
            BB = Cell_parametres::B_AXIS_Lower;
            B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        }
        else
        {
            BB = Cell_parametres::B_AXIS;
            B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        }
        double numerator = Cell_parametres::A_AXIS * BB;
        double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;

        double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;
        denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
        denominator = sqrt( denominator );
        double c_multiplicator = numerator / denominator;
        Vector3d position_updated = position * c_multiplicator;


        Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );
        Vector3d divide_position_updated = position_updated + add_to_get_segment;

        unsigned int x_index = divide_position_updated( 0 ) / this->x_width;
        unsigned int y_index = divide_position_updated( 1 ) / this->y_width;
        unsigned int z_index = divide_position_updated( 2 ) / this->z_width;

        unsigned int ID_map_index = x_index * this->neccessary_dimension * this->neccessary_dimension + y_index * this->neccessary_dimension + z_index;
        this->surface[ ID_map_index ].push_back( position_updated );

        counter++;
    }


    double difference_density = density_in_IS - density;
    Vector3d orientation_of_IS_cortical_sliding = IS_Cortical_Sl_parameter::orientation;

    if( difference_density <= 0 )
    {
	return;
    }
    else
    {
    	double area_of_IS_cor_sl = sim_of_Cell::PI * IS_Cortical_Sl_parameter::radius * IS_Cortical_Sl_parameter::radius;
    	unsigned int number_of_dynein_IS = ( unsigned int ) ( area_of_IS_cor_sl * difference_density * micro_2_to_meter_2 );
    	if( IS_Cortical_Sl_parameter::radius >= this->A_axis / 2.0 )
    	{
        	cout<<"IS_Cortical_Sl_parameter::radius >= this->A_axis / 2.0"<<endl;
        	cout<<"Surface::Surface(  double density , double width , string cortical_Sliding , double density_in_IS )"<<endl;
        	unsigned int ERROR_ID = 121632147;
        	cout<<"ERROR_ID = "<<ERROR_ID<<endl;
    	}

    	double angle_of_IS = asin( IS_Cortical_Sl_parameter::radius / this->A_axis);
    	Vector3d IS_orientation = IS_Cortical_Sl_parameter::orientation;


    	unsigned int counter_IS = 0;

    	cout<<"number_of_dynein_IS = "<<number_of_dynein_IS<<endl;

    std::uniform_real_distribution<> distribution{  ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius };
    unsigned int number_of_generator = omp_get_thread_num();

    	while( counter_IS < number_of_dynein_IS )
    	{
        	double x = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        	double y = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        	double z = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        	Vector3d position( x , y , z );

        	if( position.norm() > virtual_sphere_radius )
        	{
           		continue;
        	}


        	double BB;
        	double B_multiplicator;
        	if( z  < 0 )
        	{
            	BB = Cell_parametres::B_AXIS_Lower;
            	B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        	}
        	else
        	{
            	BB = Cell_parametres::B_AXIS;
            	B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        	}
        	double numerator = Cell_parametres::A_AXIS * BB;
        	double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;

        	double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;
        	denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
        	denominator = sqrt( denominator );
        	double c_multiplicator = numerator / denominator;
        	Vector3d position_updated = position * c_multiplicator;

        	double angle_tmp = acos( IS_orientation.dot( position_updated ) / ( position_updated.norm() * IS_orientation.norm() ) );
        	if( angle_tmp < angle_of_IS )
        	{

            		Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );
            		Vector3d divide_position_updated = position_updated + add_to_get_segment;

            		unsigned int x_index = divide_position_updated( 0 ) / this->x_width;
            		unsigned int y_index = divide_position_updated( 1 ) / this->y_width;
            		unsigned int z_index = divide_position_updated( 2 ) / this->z_width;

            		unsigned int ID_map_index = x_index * this->neccessary_dimension * this->neccessary_dimension + y_index * this->neccessary_dimension + z_index;
            		this->surface[ ID_map_index ].push_back( position_updated );

            		counter_IS++;
        	}
    	}
    }


}





Surface::Surface( string cortical_Sliding , ISCorticalSl ISCorticalSl_arg , unsigned int minus_integer )
{


    //Surface::Surface(  double density , double width , string cortical_Sliding , double density_in_IS , ISCorticalSl ISCorticalSl_arg , unsigned int minus_integer )
    this->A_axis = Cell_parametres::A_AXIS;
    this->B_axis = Cell_parametres::B_AXIS_Lower;
    this->x_width = sim_of_Cell::surface_width; //sqrt( Dynein::L_0 * Dynein::L_0 / 3.0 );
    this->y_width = sim_of_Cell::surface_width;
    this->z_width = sim_of_Cell::surface_width;
    this->density = density;
    //surface of the ellipse - so far just the surface of the sphere
    double virtual_sphere_radius = ( this->A_axis + this->B_axis ) / 2.0;
    double surface = 4.0 * sim_of_Cell::PI * virtual_sphere_radius * virtual_sphere_radius;

    double micro_2_to_meter_2 = 1e12;
    unsigned int number_of_dynein = surface * sim_of_Cell::density_of_dynein_motor_surface * micro_2_to_meter_2;
    unsigned int counter = 0;



    Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );

    std::uniform_real_distribution<> distribution{  ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius };
    unsigned int number_of_generator = omp_get_thread_num();

    while( counter < number_of_dynein )
    {
        //double x = randx( ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius );
        //double y = randx( ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius );
        //double z = randx( ( - 1.0 ) * virtual_sphere_radius ,  virtual_sphere_radius );

        double x = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        double y = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        double z = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        Vector3d position( x , y , z );

        if( position.norm() > virtual_sphere_radius )
        {
            continue;
        }


        double BB;
        double B_multiplicator;
        if( z  < 0 )
        {
            BB = Cell_parametres::B_AXIS_Lower;
            B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        }
        else
        {
            BB = Cell_parametres::B_AXIS;
            B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        }
        double numerator = Cell_parametres::A_AXIS * BB;
        double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;

        double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;
        denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
        denominator = sqrt( denominator );
        double c_multiplicator = numerator / denominator;
        Vector3d position_updated = position * c_multiplicator;


        //these two points will give me the line and the points closer to the line than cut_radius are eliminated


        Vector3d divide_position_updated = position_updated + add_to_get_segment;

        unsigned int x_index = divide_position_updated( 0 ) / this->x_width;
        unsigned int y_index = divide_position_updated( 1 ) / this->y_width;
        unsigned int z_index = divide_position_updated( 2 ) / this->z_width;

        unsigned int ID_map_index = x_index * this->neccessary_dimension * this->neccessary_dimension + y_index * this->neccessary_dimension + z_index;
        this->surface[ ID_map_index ].push_back( position_updated );

        counter++;
    }

    //filling of IS_cortical sliding
    double difference_density = sim_of_Cell::density_of_dynein_motor_Cortical_Sliding - sim_of_Cell::density_of_dynein_motor_surface;



    Vector3d orientation_of_IS_cortical_sliding = IS_Cortical_Sl_parameter::orientation;

    if( difference_density <= 0 )
    {
	return;
    }

    double area_of_IS_cor_sl = sim_of_Cell::PI * IS_Cortical_Sl_parameter::radius * IS_Cortical_Sl_parameter::radius;
    unsigned int number_of_dynein_IS = ( unsigned int ) ( area_of_IS_cor_sl * difference_density * micro_2_to_meter_2 );
    if( number_of_dynein_IS < minus_integer )
    {
	cout<<"number_of_dynein_IS < minus_integer"<<endl;
        cout<<"Surface ERROR_ID = 645158455"<<endl;
	throw("");

    }

    number_of_dynein_IS = number_of_dynein_IS - minus_integer;

    std::vector<Vector3d>  numbers_in_IS = ISCorticalSl_arg.create_IS_Cortical_Sliding_points( number_of_dynein_IS );


    for( unsigned int point_id = 0 ; point_id < numbers_in_IS.size() ; point_id ++ )
    {
	Vector3d position_updated = numbers_in_IS[ point_id ];
        Vector3d divide_position_updated = position_updated + add_to_get_segment;
        //cout<<position_updated<<endl;

        unsigned int x_index = divide_position_updated( 0 ) / this->x_width;
        unsigned int y_index = divide_position_updated( 1 ) / this->y_width;
        unsigned int z_index = divide_position_updated( 2 ) / this->z_width;

        unsigned int ID_map_index = x_index * this->neccessary_dimension * this->neccessary_dimension + y_index * this->neccessary_dimension + z_index;
        this->surface[ ID_map_index ].push_back( position_updated );

    }


}


bool Surface::dynein_surface_catch_control_tangent( Vector3d point_position , Vector3d a_point , Vector3d tangent_b , Vector3d& Catching_point )
{
    //this function control whether the map contains the point
    //if so, the points get erased and return in Catching_point

        double BB;
        double B_multiplicator;
        if( point_position( 2 ) < 0 )
        {
            BB = Cell_parametres::B_AXIS_Lower;
            B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        }
        else
        {
            BB = Cell_parametres::B_AXIS;
            B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        }
        double numerator = Cell_parametres::A_AXIS * BB;
        double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;

        double denominator = ( point_position( 0 ) * point_position( 0 ) + point_position( 1 ) * point_position( 1 ) ) * B_multiplicator;
        denominator = denominator + point_position( 2 ) * point_position( 2 ) * A_multiplicator;
        denominator = sqrt( denominator );
        double c_multiplicator = numerator / denominator;
        Vector3d position_updated = point_position * c_multiplicator;
        //I want to get the right segment
        Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );

        Vector3d position_to_get_index = position_updated + add_to_get_segment;
        unsigned int x_index = position_to_get_index( 0 ) / this->x_width;
        unsigned int y_index = position_to_get_index( 1 ) / this->y_width;
        unsigned int z_index = position_to_get_index( 2 ) / this->z_width;


        unsigned int ID_map_index = x_index * this->neccessary_dimension * this->neccessary_dimension + y_index * this->neccessary_dimension + z_index;
        try
        {
            std::vector<Vector3d> dynein_on_surface_points = this->surface.at( ID_map_index );
            for( unsigned int point_id = 0 ; point_id < dynein_on_surface_points.size() ; point_id ++ )
            {
                Vector3d dynein_position = dynein_on_surface_points.at( point_id );
                double vzdalenost = distance_point_segment( a_point , tangent_b , dynein_position );

                if( abs( vzdalenost ) < IS_Dynein_Cell_surface::dynein_on_surface_catch_radius )
                {
                    Catching_point = dynein_position;
                    this->erase_dynein_point( ID_map_index , point_id );
                    return true;
                }
            }

        }
        catch (const std::out_of_range& oor)
        {
            Catching_point = Vector3d( 0.0 , 0.0 , 0.0 );
            return false;
        }

        //this wil never happen -
        return false;
}



