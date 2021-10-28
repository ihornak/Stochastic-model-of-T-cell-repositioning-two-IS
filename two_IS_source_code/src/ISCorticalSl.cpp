#include "ISCorticalSl.h"


//Default constructor
ISCorticalSl::ISCorticalSl()
{
    
    
}

ISCorticalSl::ISCorticalSl( Vector3d axis_tmp , Vector3d plane_point_tmp , double layer_tmp )
{
    this->axis = axis_tmp;
    this->plane_point = plane_point_tmp;
    this->layer = layer_tmp; 
    //confirmation that plane_point and center of coordinates can be connected with the axis
    {
        //axis aiming at opposite direction must go from the center of coordinates to the plane_point
        Vector3d axis_tmp = this->axis * ( -1.0 );
        double ratio = this->plane_point( 0 ) /  axis_tmp( 0 );
        
        if( ( abs( this->plane_point( 1 ) - ratio * axis_tmp( 1 ) ) < 1e-10 ) && ( abs( this->plane_point( 2 ) - ratio * axis_tmp( 2 ) ) < 1e-10 )  )
        {
            cout<<" Problems in ISCorticalSl::ISCorticalSl( axis_tmp plane_point_tmp layer_tmp) "<<endl;
            cout<<" ISCorticalSl ERROR_ID = 6463416546345874"<<endl;
            throw("");
        }        
    }
    
    this->radius = sqrt( Cell_parametres::A_AXIS * Cell_parametres::A_AXIS - this->plane_point.norm() * this->plane_point.norm() );
}



//This one is used
//The cylinder of the IS cortical sliding is constructed using the centers of the front(close to IS) and the rear side
//and the radius
//The IS can have a form of a ring by setting radius_inner_tmp>  0 //So far it was not used
ISCorticalSl::ISCorticalSl( Vector3d center_front_of_IS_tmp , Vector3d center_rear_of_IS_tmp , double radius_tmp  , double radius_inner_tmp)
{
	Vector3d axis_of_IS = ( - 1.0 ) * ( center_front_of_IS_tmp - center_rear_of_IS_tmp ) / ( center_front_of_IS_tmp - center_rear_of_IS_tmp ).norm();
	this->axis = axis_of_IS;
	this->center_front_of_IS = center_front_of_IS_tmp;
	this->center_rear_of_IS = center_rear_of_IS_tmp;
	this->radius = radius_tmp;
	this->radius_inner = radius_inner_tmp;
	this->layer = ( center_front_of_IS_tmp - center_rear_of_IS_tmp ).norm();
}






ISCorticalSl::ISCorticalSl(const ISCorticalSl& other)
{
    this->axis = other.axis;
    this->center_front_of_IS = other.center_front_of_IS;
    this->center_rear_of_IS = other.center_rear_of_IS;
    this->radius = other.radius;
    this->radius_inner = other.radius_inner;
    this->layer = other.layer;

}
ISCorticalSl& ISCorticalSl::operator=( const ISCorticalSl &other )
{
    this->axis = other.axis;
    this->center_front_of_IS = other.center_front_of_IS;
    this->center_rear_of_IS = other.center_rear_of_IS;
    this->radius = other.radius;
    this->radius_inner = other.radius_inner;
    this->layer = other.layer;
    return *this;
}


//Controls whether it belongs to the IS
bool ISCorticalSl::check_ISCorticalSl_caught( Vector3d position )
{
    Vector3d point_on_IS = this->plane_point;
    Vector3d axis_of_IS = this->axis;
    double distance = distance_plane_point( axis_of_IS , point_on_IS , position );
    
    if( distance <= 0 )
    {
        return true;
    }
    else
    {
        return false;
    }
}



//Computes the distance of the segment in the IS
double ISCorticalSl::covered_distance_of_segment( Vector3d first_point , Vector3d second_point )
{

    double distance_of_segment;
    if( ( this->check_ISCorticalSl_caught( first_point ) == 0 ) && ( this->check_ISCorticalSl_caught( second_point ) == 0 ) )
    {
        return 0.0;
    }
    
    if( ( this->check_ISCorticalSl_caught( first_point ) == 1 ) && ( this->check_ISCorticalSl_caught( second_point ) == 1 ) )
    {
        distance_of_segment = ( second_point - first_point ).norm();
    }
    else if( ( this->check_ISCorticalSl_caught( first_point ) == 1 ) && ( this->check_ISCorticalSl_caught( second_point ) == 0 )  )
    {
        double distance_to_plane = abs( distance_plane_point( this->axis , this->plane_point , first_point ) );
        Vector3d orientation = second_point - first_point;
        double cos_angle = orientation.dot( this->axis ) / ( this->axis.norm() * orientation.norm() );
        distance_of_segment = distance_to_plane / cos_angle;
    }
    else 
    {
        double distance_to_plane = abs( distance_plane_point( this->axis , this->plane_point , second_point ) );
        Vector3d orientation = first_point - second_point;
        double cos_angle = orientation.dot( this->axis ) / ( this->axis.norm() * orientation.norm() );
        distance_of_segment = distance_to_plane / cos_angle;
    }
    return distance_of_segment;
}






//Returns the point of the plane
Vector3d ISCorticalSl::get_point_on_plane()
{
    return this->plane_point;
}

//Returns the axis of the IS
Vector3d ISCorticalSl::get_axis()
{
    return this->axis;
}


//Returns the center  of the front wall(close to the center of the cell)
Vector3d ISCorticalSl::get_center_front_of_IS()
{
	return this->center_front_of_IS;

}

//Returns the center  of the rear wall(distant from the center of the cell)
Vector3d ISCorticalSl::get_center_rear_of_IS()
{
	return this->center_rear_of_IS;
}

//Gets the inner radius of the IS
//In priciple, cortical sliding IS can be modelled as a ring
//This is not used: inner radius = 0
double ISCorticalSl::get_radius_inner()
{
	return this->radius_inner;
}	

//Returns the radius of the IS
double ISCorticalSl::get_radius()
{
	return this->radius;
}


//Function designed just for visualisation purposes
std::vector<Vector3d> ISCorticalSl::create_IS_Cortical_Sliding_points( unsigned int number_of_points )
{

	cout<<"Do not use this function"<<endl;
	cout<<"std::vector<Vector3d> ISCorticalSl::create_IS_Cortical_Sliding_points( unsigned int number_of_points )"<<endl;
	cout<<"This function is ok, but you are using geometrical object to create map content-> unhealthy practice"<<endl;
	cout<<"This function can be used only for the elementary visualisations"<<endl;
	throw("");
	Vector3d z_axis( 0.0 , 0.0 , 1 );
        std::vector<Vector3d> returned_points;
        double angle = acos( z_axis.dot( this->axis ) / ( this->axis.norm() * z_axis.norm() ) );
	Vector3d axis_of_rotation( 0.0 , 0.0 , 0.0  );
	if( abs( angle - sim_of_Cell::PI ) < 0.1 )
	{
		axis_of_rotation( 0 ) = 1.0; 

	}
	else
	{
		axis_of_rotation = z_axis.cross( this->axis );
        	axis_of_rotation = axis_of_rotation / axis_of_rotation.norm();
	}
	

        Quaternion<double> q;
	q = AngleAxis<double>( angle , axis_of_rotation );


	for( unsigned int number = 0 ; number < number_of_points ; number ++ )
	{

		double azimuthal_angle = rand_x( 0.0 , 2.0 * sim_of_Cell::PI );  
                double radius_tmp = triangular_distribution();
		double x = radius_tmp * cos( azimuthal_angle );
                double y = radius_tmp * sin( azimuthal_angle );

		double z = sqrt(  Cell_parametres::B_AXIS *  Cell_parametres::B_AXIS - ( x * x + y * y ) );


		Vector3d bod( x , y , z );
                bod = project_point_on_surface( bod ); 
		Vector3d rotated_point = q * bod;
		returned_points.push_back( rotated_point );
	}        

	return returned_points;
}

//Checks how many of controlled point belongs to the IS
//Whether they are inside, or in the close proximity of the IS cylinder
std::vector<Vector3d> ISCorticalSl::control_IS_Cortical_Sliding_points( std::vector<Vector3d> points_to_control , unsigned int& number_of_IS_points )
{
	double epsilon = 2e-7;
	std::vector<Vector3d> returned_points;
	unsigned int counter = 0;
	for( unsigned int point_id = 0 ; point_id < points_to_control.size() ; point_id ++  )
	{
		Vector3d bod = points_to_control[ point_id ];
		Vector3d first_line_point( 0.0 , 0.0 , 0.0 );
		Vector3d second_point = this->axis; 
                double distance = distance_point_line( first_line_point , second_point , bod );
		if (distance < this->radius + epsilon )
		{
			double cos_angle = bod.dot( this->axis ) / ( bod.norm() * this->axis.norm() ); 
			if( cos_angle > 0.0 )
			{
				counter = counter + 1;
			}			
			else
			{
				returned_points.push_back( bod );
			}
		}
		else
		{    
			returned_points.push_back( bod );
		}

	}
	number_of_IS_points = counter;
	return returned_points;
}





//Simplified control whether the point belongs to the IS
bool ISCorticalSl::control_IS_Cortical_Sliding_points_for_dynein_testing( Vector3d point )
{
	double epsilon = 1.8e-6;
	unsigned int counter = 0;
	Vector3d bod = point;
	Vector3d first_line_point( 0.0 , 0.0 , 0.0 );
	Vector3d second_point = this->axis; 
	//Gets the distance between the point and the axis
        double distance = distance_point_line( first_line_point , second_point , bod );
	if (distance < this->radius + epsilon )
	{
		double cos_angle = bod.dot( this->axis ) / ( bod.norm() * this->axis.norm() ); 
		if( cos_angle > 0.0 )
		{
			if( second_point[ 0 ] > 0 )
			{
				if( point[ 0 ] > 0 )
				{
					return true;
				}
				else
				{
					return false;				
				}
			} 
			if( second_point[ 0 ] < 0 )
			{
				if( point[ 0 ] < 0 )
				{
					return true;
				}
				else
				{
					return false;				
				}
			} 			

		}			
		else
		{
			return false;
		}
	}
	else
	{    
		return false;
	}
	return false;
}






//Checks how many of controlled point belongs to the IS
std::vector<Vector3d> ISCorticalSl::control_IS_Cortical_Sliding_points2( std::vector<Vector3d> points_to_control , unsigned int& number_in  , unsigned int& number_out )
{
	//Returned points will be projected to the cell membrane and added to the surface
	//This function works also in the scenario, where dyneins are located outside of the IS
	if( points_to_control.size() == 0 )
	{
		number_in = 0;
		number_out = 0;
		return points_to_control;
	}


	double epsilon = 6.0e-7;
	//returned jsou ty body, ktere se proste vrati a budou se projektovat na plochu
	std::vector<Vector3d> returned_points;


	unsigned int counter_1 = 0;
        unsigned int counter_2 = 0;
	for( unsigned int point_id = 0 ; point_id < points_to_control.size() ; point_id ++  )
	{
		Vector3d bod = points_to_control[ point_id ];
		Vector3d first_line_point( 0.0 , 0.0 , 0.0 );
		Vector3d second_point = this->axis; 
                double distance = distance_point_line( first_line_point , second_point , bod );
                double cos_angle = bod.dot( this->axis ) / ( bod.norm() * this->axis.norm() ); 
		if( cos_angle < 0.0 )
		{	
			returned_points.push_back( bod );
		}
		else
		{
			if ( ( distance > this->radius + epsilon ) || ( distance < this->radius ) )
			{
				//The points will be projected 	
				returned_points.push_back( bod );
			}
			//More likely belong to the IS			
			else if ( ( distance > this->radius ) && ( distance < this->radius + epsilon / 2.0 ) )
			{
				counter_1 = counter_1 + 1;
			}
			//More likely belong outside IS			
			else if ( ( distance > this->radius + epsilon / 2.0 ) && ( distance < this->radius + epsilon ) )
			{
				counter_2 = counter_2 + 1;
			}
		}
		

	}
	number_in = counter_1;
        number_out = counter_2;
	return returned_points;
}





//Function controls whether the points belong to the IS
std::vector<Vector3d> ISCorticalSl::control_IS_Cortical_Sliding_points_2_b( std::vector<Vector3d> points_to_control , std::vector<Vector3d>& inside_points )
{
	//If they do, they are saved into  inside_points
	std::vector<Vector3d> outside_points;
	if( points_to_control.size() == 0 )
	{
		return outside_points;
	}

	if( ( inside_points.size() != 0 )  )
	{
		cout<<"  ( inside_points.size() != 0 )  ) "<<endl;
		cout<<"std::vector<Vector3d> ISCorticalSl::control_IS_Cortical_Sliding_points_2_b( std::vector<Vector3d> points_to_control , std::vector<Vector3d> points_in_IS_to_project  , unsigned int& number_out )"<<endl;
		cout<<"ERROR_ID = 461646566113"<<endl;
		throw("");
	}



	//unsigned int counter = 0;
	for( unsigned int point_id = 0 ; point_id < points_to_control.size() ; point_id ++  )
	{
		Vector3d bod = points_to_control[ point_id ];
		bool answer = this->control_IS_Cortical_Sliding_points( bod );
		if( answer == false )
		{
			outside_points.push_back( bod );
		}
		if( answer == true )
		{
			inside_points.push_back( bod );
		}

	}
	return outside_points;
}




//Checks whether point belongs to the cortical sliding IS
bool ISCorticalSl::control_IS_Cortical_Sliding_points( Vector3d point )
{
	double epsilon = 1.0e-6;
	unsigned int counter = 0;
	Vector3d bod = point;
	Vector3d first_line_point( 0.0 , 0.0 , 0.0 );
	Vector3d second_point = this->axis; 
        double distance = distance_point_line( first_line_point , second_point , bod );
	if (distance < this->radius + epsilon )
	{
		double cos_angle = bod.dot( this->axis ) / ( bod.norm() * this->axis.norm() ); 
		if( cos_angle > 0.0 )
		{
			return true;
		}			
		else
		{
			return false;
		}
	}
	else
	{    
		return false;
	}
	return false;
}


//Controls whether points belong to the IS
std::vector<Vector3d> ISCorticalSl::control_IS_Cortical_Sliding_points_2_b_for_dynein_testing( std::vector<Vector3d> points_to_control , std::vector<Vector3d>& inside_points )
{
	//If they do, they go to inside_points
	//outside_points otherwise
	std::vector<Vector3d> outside_points;
	if( points_to_control.size() == 0 )
	{
		return outside_points;
	}

	if( ( inside_points.size() != 0 )  )
	{
		cout<<"  ( inside_points.size() != 0 )  ) "<<endl;
		cout<<"std::vector<Vector3d> ISCorticalSl::control_IS_Cortical_Sliding_points_2_b( std::vector<Vector3d> points_to_control , std::vector<Vector3d> points_in_IS_to_project  , unsigned int& 			number_out )"<<endl;
		cout<<"ERROR_ID = 461646566113"<<endl;
		throw("");
	}

	for( unsigned int point_id = 0 ; point_id < points_to_control.size() ; point_id ++  )
	{
		Vector3d bod = points_to_control[ point_id ];
		bool answer = this->control_IS_Cortical_Sliding_points_for_dynein_testing( bod );
		if( answer == false )
		{
			outside_points.push_back( bod );
		}
		if( answer == true )
		{
			inside_points.push_back( bod );
		}

	}
	return outside_points;
}






//Controls the distance from the axis
bool ISCorticalSl::control_IS_Cortical_Sliding_points_inside_meze( Vector3d point )
{
        double epsilon = 6.0e-7;;
	unsigned int counter = 0;
	Vector3d bod = point;
	Vector3d first_line_point( 0.0 , 0.0 , 0.0 );
	Vector3d second_point = this->axis; 
	double cos_angle = bod.dot( this->axis ) / ( bod.norm() * this->axis.norm() ); 
	
	if( cos_angle < 0.0 )
	{
			cout<<"else"<<endl;
			cout<<"bool ISCorticalSl::control_IS_Cortical_Sliding_points( Vector3d point )"<<endl;
			cout<<"ERROR_ID = 3131353"<<endl;
			throw("");
	}
        double distance = distance_point_line( first_line_point , second_point , bod );	
	if( distance <=  this->radius )
	{
		return true;
	}
	else
	{
		return false;
	}
}








//Projects the point on the cell membrane
Vector3d ISCorticalSl::project_point_on_surface( Vector3d position )
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









//Controls whether the points are in the IS or at the margins
std::vector<Vector3d> ISCorticalSl::control_of_inside_and_margins( std::vector<Vector3d> points_to_control , std::vector<Vector3d>& points_margins )
{
        //Control of the input
	std::vector<Vector3d> inside_points;
	if( ( points_to_control.size() == 0 )  )
	{
		return inside_points;
	}
	if( ( points_margins.size() != 0 )  )
	{
		cout<<"  ( points_margins.size() != 0 )  ) "<<endl;
		cout<<"std::vector<Vector3d> ISCorticalSl::control_of_inside_and_margins( std::vector<Vector3d> points_to_control , std::vector<Vector3d> points_in_IS_to_project  , unsigned int& number_out )"<<endl;
		cout<<"ERROR_ID = 6846468"<<endl;
		throw("");
	}


	for( unsigned int point_id = 0 ; point_id < points_to_control.size() ; point_id ++  )
	{
		Vector3d bod = points_to_control[ point_id ];
    		bool answer = this->control_IS_Cortical_Sliding_points_inside_meze( bod );
		if( answer == true )
		{
			inside_points.push_back( bod );
		}
		else
		{
			points_margins.push_back( bod );
		}

	}
	return inside_points;
}




ISCorticalSl::~ISCorticalSl()
{
        
}
