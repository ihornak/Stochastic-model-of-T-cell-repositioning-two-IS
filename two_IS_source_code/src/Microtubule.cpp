/*
 * Microtubule.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: hornak
 */

#include "Microtubule.h"

Microtubule::Microtubule()
{
	this->cell_id = 0;
	//simple default constructor - direct microtubule
	// all beads on x axis
	this->microtubule_id = 10000;
	this->polygon = 10000;
	this->restDistancePoints = sim_of_Cell::resting_distance;
        this->MTOC_point = 1;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	//radius has no meaning - the cell do not restrict Microtubule
	this->dydein_index = 0.0;
	this->effective_friction = 1e-7;

	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 0.0 , 0.0 , 0.0 );


	if( this->numberOfPoints < 1  )
	{
		cout<<" this->numberOfPoints < 1 in Microtubule::Microtubule() "<<endl;
		throw("");
	}

	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		this->coordinates( 3 * i , 0 ) = i * this->restDistancePoints;
		this->coordinates( 3 * i + 1 , 0 ) = 0.0;
		this->coordinates( 3 * i + 2 , 0 ) = 0.0;
	}

	//this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
    this->set_growth_index( 0 );


    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	
    	
}



Microtubule::Microtubule( unsigned int ID )
{
	this->cell_id = 0;
	//simple default constructor - direct microtubule just with id
	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
        this->MTOC_point = 1;
	this->polygon = 10000;
	//radius has no meaning - the cell do not restrict Microtubule
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
	this->effective_friction = this->calculateEffectiveFriction_Howard();
	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 0.0 , 0.0 , 0.0 );
	if( this->numberOfPoints < 1  )
	{
		cout<<" this->numberOfPoints < 1 in Microtubule::Microtubule( unsigned int ID ) "<<endl;
		throw("");
	}

	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );

	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		this->coordinates( 3 * i , 0 ) = i * this->restDistancePoints;
		this->coordinates( 3 * i + 1 , 0 ) = 0.0;
		this->coordinates( 3 * i + 2 , 0 ) = 0.0;
	}
	//this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
    this->set_growth_index( 0 );
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	
}





Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points  )
{
    	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points
        //First segment can be shorter, since it is part of the MTOC

	//This can be potentially used as the cell(process) identificator
	this->cell_id = 0;
    	if( MTOC_point_arg == 0 )
    	{
        	cout<<"MTOC_point_arg == 0"<<endl;
        	cout<<"ERROR_ID = 689645684316464"<<endl;
        	throw("");
    	}
    	if( second_MTOC_point_arg == 0 )
    	{
        	cout<<"MTOC_point_arg == 0"<<endl;
        	cout<<"ERROR_ID = 8987625148688"<<endl;
        	throw("");
    	}

	//Resting distance of the segments
    	this->restDistancePoints = sim_of_Cell::resting_distance;
    	//bending rigidity of discretized case
    	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
    	this->numberOfPoints = number_of_points;
    	
	//These are the identificators of the uni
	//Unused in the current version of the code where the MTOC will be 3d structure    	
    	this->microtubule_id = ID;
    	this->polygon = poly;
    	this->side = side_arg;
    	this->MTOC_point = MTOC_point_arg;
        //Indexes of the MTOC points, to which the microtubule is attached    	
    	this->second_MTOC_point = second_MTOC_point_arg;
    	//Dynein index determines whether the microtubule is unattached, attached to cortical sliding or capture-shrinkage dynein
    	//0 - microtubule is unattached
	this->dydein_index = 0.0;

	this->effective_friction = this->calculateEffectiveFriction_Howard();
	//the length of the segment in the MTOC
	this->restDistancePoints_first = ( first_Point - second_Point ).norm();

    	if( abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10  )
    	{
        	cout<<"( first_Point - second_Point ).norm() = "<<( first_Point - second_Point ).norm()<<endl;
        	cout<<"this->restDistancePoints_first = "<<this->restDistancePoints_first<<endl;
        	cout<<"abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10"<<endl;
        	cout<<"ERROR_ID = 61543554435"<<endl;
        	throw("");
    	}
	//Constrols bad input
	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );
    	double a_axis = Cell_parametres::A_AXIS;
    	double b_axis = Cell_parametres::B_AXIS;
    	//confirm that first bead lies in the ellipsoid - if not throw error
	bool bool_index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );

	if( bool_index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"bool_index = "<<bool_index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=46461343154127"<<endl;
		throw("");
	}


	//this is direction of microtubule in x - y plane - does not change for entire microtubule
	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	//axis that I will use to rotate the vectors
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
	{
		this->coordinates( 3 + dimension , 0 ) = second_Point( dimension );
	}


	//Last point and orientation will be changed as we go from
	Vector3d last_point = second_Point;
	Vector3d last_orientation = orientation;
	//The microtubule is created
	//Bacially, the orientation of the last segment is rotated by one angle to create orientation of the new segment
	//If the new bead is outside the boundary of the cell, the segment is rotated once more
	//It assures that the microtubule copy the cell membrane	
	for( unsigned int index = 2 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		//The coordinates are writen into matrix
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{
			//I aplly quaternion to rotate
			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//Rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;

			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    //Lengths of the segments of microtubules are the same
    //The exception is the first segment, since it is the part of the MTOC and is fitted to the two points of the MTOC	
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    //Bending matrix is unused - for the future simulation
    this->set_bending_matrix( );
    //Grow index determines whether microtubule is growing, shrinking....
    //At the beginning, it is stable  - 0
    this->set_growth_index( 0 );
    //Sets the length of the microtubule outside the MTOC
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    //Sets the lengths of segments between the beads of the microtubule
    this->set_lenght_of_tangents_micro_0();
}















Microtubule::Microtubule( const Microtubule & tmp )
{

	this->restDistancePoints = tmp.getRestDist();
	this->kappa =  tmp.getKappa();// 2.2e-25 v pohode  2.2e-26
	this->numberOfPoints = tmp.getNumberOfPoints();
	this->microtubule_id = tmp.microtubule_id;
	this->effective_friction = tmp.effective_friction;
	this->IS_position_catching = tmp.IS_position_catching;

	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates = tmp.getCoordinates();
	this->polygon = tmp.polygon;
	this->dydein_index = tmp.dydein_index;


    //this->plane_axis_point = tmp.plane_axis_point;
    //this->Dynein_on_surface = tmp.Dynein_on_surface;
    this->Dynein_motors_2 = tmp.Dynein_motors_2;
    this->MTOC_point = tmp.MTOC_point;
    this->second_MTOC_point = tmp.second_MTOC_point;
    this->restDistancePoints_first = tmp.restDistancePoints_first;
    this->lenght_of_tangents = tmp.lenght_of_tangents;
    this->b_MATRIX = tmp.b_MATRIX;
    this->Kronecker_Delta = tmp.Kronecker_Delta;   
    this->growth_index = tmp.growth_index;
    this->lenght_of_micro_0_outside_MTOC = tmp.lenght_of_micro_0_outside_MTOC;	    	
    this->lenght_of_tangents_micro_0 = tmp.lenght_of_tangents_micro_0;
    this->cell_id = tmp.cell_id;	
}



//Sets the coordinates of the microtubule
void Microtubule::setCoordinates( MatrixXd &coordinatesTmp )
{
	//Controls mistakes in the input
	if( ( coordinatesTmp.rows() !=  3 * this->numberOfPoints ) || ( ( coordinatesTmp.cols() != 1 ) ) )
	{
		cout<<"( ( coordinatesTmp.rows() !=  3 * this->numberOfPoints ) || ( ( coordinatesTmp.cols() != 1 ) ) in Microtubule::setCoordinates( MatDoub setCoordinates )"<<endl;
		cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		throw("");
	}
	this->coordinates = coordinatesTmp;
}


//Sets the position of one bead
void Microtubule::setPoint( unsigned int number_of_point ,  Vector3d position )
{

    if( number_of_point >= this->getNumberOfPoints() )
    {
        cout<<"( number_of_point >= this->getNumberOfPoints() )"<<endl;
        cout<<"Microtubule::setPoint( unsigned int number_of_point ,  Vector3d position )"<<endl;
        throw("");
    }

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        //cout<<"position( dimension ) = "<<position( dimension )<<endl;
        this->coordinates( 3 * number_of_point + dimension , 0 ) = position( dimension );
    }
}






//Returns the coordinates of the microtubule
MatrixXd Microtubule::getCoordinates()  const
{
	return this->coordinates;
}

//Returns the length of the segment
double Microtubule::getRestDist() const
{
	return this->restDistancePoints;
}

//Returns the bending rigidity of the discretized polymer
double Microtubule::getKappa()  const
{
	return this->kappa;
}








//Returns the position, where the microtubule depolimerizes
Vector3d Microtubule::get_IS_position_catching( ) const
{
	return this->IS_position_catching;
}

//Returns the effective friction of the microtubule bead
double Microtubule::get_effective_friction()
{
	return this->effective_friction;
}

double Microtubule::get_effective_friction_whole_microtubule()
{
	return this->effective_friction * ( (double) this->numberOfPoints );
}



//Returns the identifier of the microtubule
unsigned int Microtubule::getID() const
{
	return this->microtubule_id;
}

//Returns the side of the MTOC where the microtubule is anchored
unsigned int Microtubule::getSide()  const
{
    return this->side;
}


//Returns the number of MTOC polygon, to which the microtubule is attached 
unsigned int Microtubule::get_polygon_number() const
{
	return this->polygon;
}


//Returns the MTOC point from which the microtubule is sprouting 
unsigned int Microtubule::get_MTOC_point()
{
    return this->MTOC_point;
}

//Return rear points of the MTOC
unsigned int Microtubule::get_MTOC_opposite_point()
{
    return this->second_MTOC_point;
}


//Determines whether the microtubule is free, or attach to capture-shrinkage or cortical sliding dyneins
unsigned int Microtubule::get_dynein_index() const
{
	return this->dydein_index;
}


//Sets the dynein index
void Microtubule::set_dynein_index( unsigned int new_value )
{
	this->dydein_index = new_value;
}

//Sets the point, where the microtubule is depolimerized
void Microtubule::set_IS_position_catching( Vector3d IS_vector )
{
	this->IS_position_catching = IS_vector;
}


 //Overloading of the operator
Microtubule& Microtubule::operator=( const Microtubule &tmp )
{
	if( this == &tmp )
	{
		cout<<" this == &tmp  in Microtubule::operator=( const Microtubule &tmp )"<<endl;
		throw("");
	}

	this->restDistancePoints = tmp.getRestDist();
	this->kappa =  tmp.getKappa();// 2.2e-25 v pohode  2.2e-26
	this->numberOfPoints = tmp.getNumberOfPoints();
	this->coordinates = tmp.getCoordinates();
	this->dydein_index = tmp.get_dynein_index();
	this->IS_position_catching = tmp.get_IS_position_catching();
	this->effective_friction = tmp.effective_friction;

    this->microtubule_id = tmp.microtubule_id;
    this->side = tmp.side;
    this->polygon = tmp.get_polygon_number();
	this->MTOC_point = tmp.MTOC_point;
    this->second_MTOC_point = tmp.second_MTOC_point;
    //this->plane_axis_point = tmp.plane_axis_point;
    //this->Dynein_on_surface = tmp.Dynein_on_surface;
    this->restDistancePoints_first = tmp.restDistancePoints_first;
    this->lenght_of_tangents = tmp.lenght_of_tangents;
    this->b_MATRIX = tmp.b_MATRIX;
    this->Kronecker_Delta = tmp.Kronecker_Delta;
    this->growth_index = tmp.growth_index;
    this->lenght_of_micro_0_outside_MTOC = tmp.lenght_of_micro_0_outside_MTOC;
    this->lenght_of_tangents_micro_0 = tmp.lenght_of_tangents_micro_0; 	
    this->cell_id = tmp.cell_id;
    return *this;
}


//Get number of mirotubule beads
unsigned int Microtubule::getNumberOfPoints()  const
{
	return this->numberOfPoints;
}


//Calculate effective friction  using mentioned publication  David Leith - Aerosol Science and Technology
double Microtubule::calculateEffectiveFriction(  )
{
	//radius of the cell is always 1.25e-8
	//Stokes law is used
	//FrictionForce = 3 * PI * viscosity * diametr * K
	//K - dynamic shape factor

	//Drag on Nonspherical Objects - David Leith - Aerosol Science and Technology

	double lenghtOfSegment = this->getRestDist();
	double viscosity =  sim_of_Cell::viscosity;  // 9.3 * 10.0

    double radius = 1.25e-8;
    double tmp = 3.0 / 4.0 * radius * radius * lenghtOfSegment;
    double d_v = 2.0 * std::pow( tmp , 1.0 / 3.0);
    double d_u = 2.0 * sqrt( 2.0 * radius * lenghtOfSegment / 3.14159265 );
    double d_s = sqrt( 2.0 * radius * radius + radius * lenghtOfSegment );			//tady dvojka neni umyslne
    double K = 1.0 / 3.0 * d_u / d_v + 2.0 / 3.0 * d_s / d_v;

    double effectiveFriction = 3.0 * sim_of_Cell::PI * sim_of_Cell::multiply_friction_constant * viscosity * d_v * K;
    return effectiveFriction;
}

//Caculates effective friction Howard, J., 2001. Mechanics of Motor Proteins and the Cytoskeleton
double Microtubule::calculateEffectiveFriction_Howard(  )
{

    double lenghtOfSegment = this->getRestDist();
    double viscosity =  sim_of_Cell::viscosity;  // 9.3 * 10.0
    double radius = 1.25e-8;

    double nominator = 4.0 * sim_of_Cell::PI * viscosity * lenghtOfSegment;
    double denominator = log( lenghtOfSegment / ( 2.0 * radius ) ) + 0.84;

    double effectiveFriction = nominator / denominator;

    return effectiveFriction;
}




//Get the point determined by the index
Vector3d Microtubule::getPoint( unsigned int index )  const
{
	//Checks the input
	if( index >= this->numberOfPoints  )
	{

		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		cout<<" this->getNumberOfPoints() = "<<this->getNumberOfPoints()<<endl;
		cout<<"index = "<<index<<endl;
		cout<<"index >= this->numberOfPoints = in Vector3d Microtubule::getPoint( unsigned int index ) "<<endl;
		throw("");
	}


	Vector3d result( 0.0 , 0.0 , 0.0 );
	for( unsigned int i = 0 ; i < 3 ; i ++ )
	{
		result( i ) = this->coordinates( 3 * ( index ) + i , 0 );
	}

	return result;
}

//Returns the last point of the microtubule
Vector3d Microtubule::get_last_Point( )  const
{
    return this->getPoint( this->getNumberOfPoints() - 1 );
}





//Returns the segment connecting two beads
Vector3d Microtubule::getTangent2( unsigned int index )  const
{
	//Controls the input
	//If the input is wrong, the simulation is stopped
	if( index >= this->numberOfPoints - 1 )
	{

        	cout<<"this->dydein_index = "<<this->dydein_index<<endl;
        	cout<<" number of points in microtubule = "<<this->getNumberOfPoints()<<endl;
        	cout<<"index = "<<index<<endl;
		cout<<"index >= this->numberOfPoints - 1 in Microtubule::getTangent2( unsigned int index )"<<endl;
        	int ERROR_ID = 1411;
		cout<<" Microtubule ERROR_ID = "<<ERROR_ID<<endl;
        	double lenght = this->get_lenght_of_microtubule();
        	cout<<"............................................."<<endl;
        	cout<<"lenght = "<<lenght<<endl;
        	if( ( this->dydein_index == 9 ) || ( this->dydein_index == 20 ) )
        	{
            		cout<<"condition"<<endl;
            		cout<<"this->Dynein_motors_2.size() = "<<this->Dynein_motors_2.size()<<endl;
            		for( unsigned int i = 0 ; i < this->Dynein_motors_2.size() ; i ++ )
            		{
                		std::pair < Vector3d , double > pair_tmp = Dynein_motors_2.at( i );
                		double abscissa = std::get<1>( pair_tmp );
                		if( abscissa > lenght )
                		{
                    			cout<<"abscissa = "<<abscissa<<endl;
                		}
           		}
        	}

        	throw ERROR_ID;

	}

	Vector3d tangent( 0.0 , 0.0 , 0.0 );
	for( unsigned int j = 0 ; j < 3 ; j ++ )
	{
		tangent( j ) = this->coordinates( 3 * ( index + 1 ) + j , 0 ) - this->coordinates( 3 * index + j , 0 );
	}
	return tangent;
}


//Gets the tangent connecting the last two points of the microtubule
Vector3d Microtubule::get_last_Tangent( )  const
{
	return getTangent2( this->getNumberOfPoints() - 2 );
}

//Prolongs the last segment of the microtubule by additional_lenght
void Microtubule::prolong_Last_tangent(  double additional_lenght )
{
	double first =  this->get_last_Tangent( ).norm();
	//Controls the input
	if( additional_lenght < 0 )
	{
		cout<<"additional_lenght < 0 "<<endl;
		cout<<"additional_lenght =  "<<additional_lenght<<endl;
		cout<<"Vector3d Microtubule::prolong_Last_tangent(  double additional_lenght ) "<<endl;
		cout<<"ERROR_ID = 3156131561361"<<endl;
		throw(" ");
	} 
	if( additional_lenght > this->getRestDist() )
	{
		cout<<" additional_lenght > this->getRestDist()  "<<endl;
		cout<<"additional_lenght =  "<<additional_lenght<<endl;
		cout<<"Vector3d Microtubule::prolong_Last_tangent(  double additional_lenght )"<<endl;
		cout<<"ERROR_ID = 130313215661"<<endl;
		throw(" ");
	}
	Vector3d previous_to_last = this->getPoint( this->getNumberOfPoints() - 2 );
	Vector3d last_tangent = this->get_last_Tangent( );
	Vector3d new_last_tangent = ( last_tangent / last_tangent.norm() ) * ( last_tangent.norm() + additional_lenght );
	Vector3d new_last_point = previous_to_last + new_last_tangent;

	this->setPoint(this->getNumberOfPoints() - 1 ,  new_last_point );
    	this->set_lenght_of_tangents("prolong_Last_tangent");


}


//Shortens the last segment of the microtubule
void Microtubule::shorten_Last_tangent(  double shortening )
{

	double first =  this->get_last_Tangent( ).norm();
	//Controls the input
	if( shortening > this->getRestDist() * 4.0 / 5.0 )
	{
		cout<<"shortening < this->getRestDist() / 5.0"<<endl;
		cout<<"shortening =  "<<shortening<<endl;
		cout<<"Vector3d Microtubule::shorten_Last_tangent(  double additional_lenght ) "<<endl;
		cout<<"ERROR_ID = 9846514566"<<endl;
		throw(" ");
	} 
	if( shortening > this->getRestDist() )
	{
		cout<<" shortening > this->getRestDist()  "<<endl;
		cout<<"shortening =  "<<shortening<<endl;
		cout<<"Vector3d Microtubule::shorten_Last_tangent(  double additional_lenght )"<<endl;
		cout<<"ERROR_ID = 949863513135646841"<<endl;
		throw(" ");
	}
	if( this->getNumberOfPoints() != 3 )
	{
		cout<<" this->getNumberOfPoints() != 3 "<<endl;
		cout<<" void Microtubule::shorten_Last_tangent(  double shortening ) "<<endl;
		cout<<"ERROR_ID = 146124631543613514"<<endl;
		throw("");
	}


	Vector3d previous_to_last = this->getPoint( 1 );
	Vector3d last_tangent = this->getTangent2( 1 );
	Vector3d new_last_tangent = ( last_tangent / last_tangent.norm() ) * ( last_tangent.norm() - shortening );
	if( new_last_tangent.norm() < dynamic_instability::coefficient_4 )
	{
		new_last_tangent = new_last_tangent / new_last_tangent.norm() * ( new_last_tangent.norm() + dynamic_instability::polymerization_constant * 0.5 ); 
	}
	Vector3d new_last_point = previous_to_last + new_last_tangent;
	this->setPoint( 2 ,  new_last_point );



    	this->set_lenght_of_tangents("prolong_Last_tangent");
	if( this->get_last_Tangent( ).norm() < dynamic_instability::coefficient_4 )
	{
		cout<<"if( this->get_last_Tangent( ) < dynamic_instability::coefficient_4 )"<<endl;
		cout<<endl;
	}



	//cout<<"this->get_last_Tangent( ) = "<<this->get_last_Tangent( ).norm() - ( first + additional_lenght )  <<endl;

}












//Gets the length of the microtubule
double Microtubule::get_lenght_of_microtubule() const
{
    double lenght_of_microtubule = 0;
    for( unsigned int segment_id = 0 ; segment_id < this->getNumberOfPoints() - 1 ; segment_id ++ )
    {
        Vector3d tangent;
        try
        {
            tangent = this->getTangent2( segment_id );
        }
        catch( int e )
        {
            cout<<"double Microtubule::get_lenght_of_microtubule() const"<<endl;
        }
        //cout<<tangent.norm()<<endl;
        lenght_of_microtubule = lenght_of_microtubule + tangent.norm();
    }
    return lenght_of_microtubule;
}

//Gets the length of the microtubule outside of MTOC
double Microtubule::get_lenght_of_microtubule_outside_MTOC() const
{
    //The first segment is the part of the MTOC
    if( this->getNumberOfPoints() <= 2 )
    {
	return 0;
    }
    double lenght_of_microtubule = 0;
    //The first segment is the part of the MTOC
    for( unsigned int segment_id = 1 ; segment_id < this->getNumberOfPoints() - 1 ; segment_id ++ )
    {
        Vector3d tangent;
        try
        {
            tangent = this->getTangent2( segment_id );
        }
        catch( int e )
        {
            cout<<"double Microtubule::get_lenght_of_microtubule() const"<<endl;
        }
        lenght_of_microtubule = lenght_of_microtubule + tangent.norm();
    }
    return lenght_of_microtubule;
}




double Microtubule::get_distance_to_lower_bead_with_index( unsigned int index ) const
{
    if( index < this->numberOfPoints - 1 )
    {
        double lenght = 0;
        for( unsigned int counter = 0 ; counter < index ; counter ++ )
        {
            lenght = lenght + this->getTangent2( counter ).norm();
        }
        return lenght;
    }
    else if( index == this->numberOfPoints - 1 )
    {
        return this->get_lenght_of_microtubule();
    }
    else
    {
        cout<<"index > this->numberOfPoints - 1 "<<endl;
        cout<<"double Microtubule::get_distance_to_lower_bead_with_index( unsigned int index ) const"<<endl;
        unsigned int ERROR_ID = 981468;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }



}


void Microtubule::set_lenght_after_catching( double lenght_after_catching_arg )
{
    this->lenght_after_catching = lenght_after_catching_arg;
}

double Microtubule::get_lenght_after_catching(  )
{
    return this->lenght_after_catching;
}


//Creates the matrix expresing derivatives of constraints
void Microtubule::getN_i_mu( MatrixXd &n_i_mu )
{
	//Brownian Dynamics algorithm for bead-rod semiflexible chain with anisotropic friction
	//Montesi Morse Pasquali
	for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		Vector3d tangent = this->getTangent2( i );
		tangent = tangent * ( 1.0 / tangent.norm() );
		for( unsigned int j = 0 ; j < 6 ; j ++ )
		{
			int tmp = 3 * i;
			if( j < 3)
			{
				n_i_mu( tmp + j , i ) = - tangent( j );
			}
			else
			{
				n_i_mu( tmp + j , i ) =  tangent( j - 3 );
			}
		}
	}
}











//Creates the matrixes for the projection of forces
void Microtubule::getMatrixes( MatrixXd &inv_Matrix , MatrixXd &projection_Matrix )
{
	//Brownian Dynamics algorithm for bead-rod semiflexible chain with anisotropic friction
	//Montesi Morse Pasquali
	//Projection matrix is expressed in equation 16
	//projectionMatrix has to be Kronecker Delta			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	MatrixXd scalars =  MatrixXd::Zero(this->numberOfPoints - 2 , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints - 2 ; i ++ )
	{
		scalars( i , 0 ) = this->getTangent2( i ).dot(this->getTangent2( i + 1 ) ) / ( this->getTangent2( i ).norm() * this->getTangent2( i + 1 ).norm() );
	}

	unsigned int coordinates = 3.0 * this->numberOfPoints;
	unsigned int bonds = this->numberOfPoints - 1;
        MatrixXd n_i_mu2 = MatrixXd::Zero( coordinates , bonds );
        getN_i_mu( n_i_mu2 );
        MatrixXd n_i_mu2Transpose = n_i_mu2.transpose();


	MatrixXd G_uv2 = MatrixXd::Zero( bonds , bonds );
	for( unsigned int index = 0 ; index < bonds ; index ++ )
	{
		G_uv2( index , index ) = 2;
		if( index != 0 )
		{
			G_uv2( index - 1 , index ) = - scalars( index - 1 , 0 );
			G_uv2( index , index - 1 ) = G_uv2( index - 1 , index );
		}
	}

	inv_Matrix = G_uv2.inverse();
	MatrixXd tmp = inv_Matrix * n_i_mu2Transpose;
	tmp = n_i_mu2 * tmp;


	projection_Matrix = projection_Matrix - tmp;
}




















void Microtubule::set_bending_matrix( )
{

	MatrixXd MATRIX = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );


	if( this->numberOfPoints == 3 )
	{
		//cout<<"AAAAAAAAAAAAAAAAAAA"<<endl;
		MATRIX( 0 , 0 ) = -1.0;
		MATRIX( 0 , 3 ) = 2.0;
                MATRIX( 0 , 6 ) = -1.0;

		MATRIX( 1 , 1 ) = -1.0;
		MATRIX( 1 , 4 ) = 2.0;
                MATRIX( 1 , 7 ) = -1.0;


		MATRIX( 2 , 2 ) = -1.0;
		MATRIX( 2 , 5 ) = 2.0;
                MATRIX( 2 , 8 ) = -1.0;

		//second bead
		MATRIX( 3 , 0 ) = 2.0;
		MATRIX( 3 , 3 ) = -4.0;
                MATRIX( 3 , 6 ) = 2.0;

		MATRIX( 4 , 1 ) = 2.0;
		MATRIX( 4 , 4 ) = -4.0;
                MATRIX( 4 , 7 ) = 2.0;

		MATRIX( 5 , 2 ) = 2.0;
		MATRIX( 5 , 5 ) = -4.0;
                MATRIX( 5 , 8 ) = 2.0;

		//third bead
		MATRIX( 6 , 0 ) = -1.0;
		MATRIX( 6, 3 ) = 2.0;
                MATRIX( 6 , 6 ) = -1.0;

		MATRIX( 7 , 1 ) = -1.0;
		MATRIX( 7, 4 ) = 2.0;
                MATRIX( 7 , 7 ) = -1.0;

		MATRIX( 8 , 2 ) = -1.0;
		MATRIX( 8, 5 ) = 2.0;
                MATRIX( 8 , 8 ) = -1.0;

		//this->b_MATRIX = MATRIX;


	}
	else if( this->numberOfPoints == 4 )
	{
		//MatrixXd bending_matrix_4 = MatrixXd::Zero( 12 , 12 );

		//first bead
		MATRIX( 0 , 0 ) = -1.0;
		MATRIX( 0 , 3 ) = 2.0;
                MATRIX( 0 , 6 ) = -1.0;


		MATRIX( 1 , 1 ) = -1.0;
		MATRIX( 1 , 4 ) = 2.0;
                MATRIX( 1 , 7 ) = -1.0;


		MATRIX( 2 , 2 ) = -1.0;
		MATRIX( 2 , 5 ) = 2.0;
                MATRIX( 2 , 8 ) = -1.0;

		//second bead
		MATRIX( 3 , 0 ) = 2.0;
		MATRIX( 3 , 3 ) = -5.0;
                MATRIX( 3 , 6 ) = 4.0;
                MATRIX( 3 , 9 ) = -1.0;

		MATRIX( 4 , 1 ) = 2.0;
		MATRIX( 4 , 4 ) = -5.0;
                MATRIX( 4 , 7 ) = 4.0;
                MATRIX( 4 , 10 ) = -1.0;


		MATRIX( 5 , 2 ) = 2.0;
		MATRIX( 5 , 5 ) = -5.0;
                MATRIX( 5 , 8 ) = 4.0;
                MATRIX( 5 , 11 ) = -1.0;


		//third bead
		MATRIX( 6 , 0 ) = -1.0;
		MATRIX( 6 , 3 ) = 4.0;
                MATRIX( 6 , 6 ) = -5.0;
                MATRIX( 6 , 9 ) = 2.0;

		MATRIX( 7 , 1 ) = -1.0;
		MATRIX( 7 , 4 ) = 4.0;
                MATRIX( 7 , 7 ) = -5.0;
                MATRIX( 7 , 10 ) = 2.0;


		MATRIX( 8 , 2 ) = -1.0;
		MATRIX( 8 , 5 ) = 4.0;
                MATRIX( 8 , 8 ) = -5.0;
                MATRIX( 8 , 11 ) = 2.0;



		//fourth bead
		MATRIX( 9 , 3 ) = -1.0;
		MATRIX( 9 , 6 ) = 2.0;
                MATRIX( 9 , 9 ) = -1.0;

		MATRIX( 10 , 4 ) = -1.0;
		MATRIX( 10 , 7 ) = 2.0;
                MATRIX( 10 , 10 ) = -1.0;

		MATRIX( 11 , 5 ) = -1.0;
		MATRIX( 11 , 8 ) = 2.0;
                MATRIX( 11 , 11 ) = -1.0;
		//this->b_MATRIX = MATRIX;


	}
	else if( this->numberOfPoints >= 5 )
	{
		unsigned int N_beads = this->numberOfPoints;
//		MatrixXd MATRIX = MatrixXd::Zero( 3 * N_beads , 3 * N_beads );
		MATRIX( 0 , 0 ) = -1.0;
		MATRIX( 0 , 3 ) = 2.0;
                MATRIX( 0 , 6 ) = -1.0;


		MATRIX( 1 , 1 ) = -1.0;
		MATRIX( 1 , 4 ) = 2.0;
                MATRIX( 1 , 7 ) = -1.0;


		MATRIX( 2 , 2 ) = -1.0;
		MATRIX( 2 , 5 ) = 2.0;
                MATRIX( 2 , 8 ) = -1.0;

		//second bead
		MATRIX( 3 , 0 ) = 2.0;
		MATRIX( 3 , 3 ) = -5.0;
                MATRIX( 3 , 6 ) = 4.0;
                MATRIX( 3 , 9 ) = -1.0;

		MATRIX( 4 , 1 ) = 2.0;
		MATRIX( 4 , 4 ) = -5.0;
                MATRIX( 4 , 7 ) = 4.0;
                MATRIX( 4 , 10 ) = -1.0;


		MATRIX( 5 , 2 ) = 2.0;
		MATRIX( 5 , 5 ) = -5.0;
                MATRIX( 5 , 8 ) = 4.0;
                MATRIX( 5 , 11 ) = -1.0;


		//last bead
		MATRIX( 3 * ( N_beads - 1 ) , 3 * ( N_beads - 3 ) ) = -1.0;
		MATRIX( 3 * ( N_beads - 1 ) , 3 * ( N_beads - 2 )  ) = 2.0;
		MATRIX( 3 * ( N_beads - 1 ) , 3 * ( N_beads - 1 )  ) = -1.0;

		MATRIX( 3 * ( N_beads - 1 ) + 1 , 3 * ( N_beads - 3 ) + 1 ) = -1.0;
		MATRIX( 3 * ( N_beads - 1 ) + 1 , 3 * ( N_beads - 2 ) + 1 ) = 2.0;
		MATRIX( 3 * ( N_beads - 1 ) + 1 , 3 * ( N_beads - 1 ) + 1 ) = -1.0;

		MATRIX( 3 * ( N_beads - 1 ) + 2 , 3 * ( N_beads - 3 ) + 2 ) = -1.0;
		MATRIX( 3 * ( N_beads - 1 ) + 2 , 3 * ( N_beads - 2 ) + 2 ) = 2.0;
		MATRIX( 3 * ( N_beads - 1 ) + 2 , 3 * ( N_beads - 1 ) + 2 ) = -1.0;

		for( unsigned int bead_id = 0 ; bead_id < N_beads ; bead_id ++ )
		{
			//cout<<"bead_id = "<<bead_id<<endl;

			if( ( bead_id > 1 ) && ( bead_id < N_beads - 2 ) )
			{
				//cout<<"AAAAAAAAAAAAAAAAAAAAAa = "<<bead_id<<endl;

				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id - 2 ) ) = - 1.0;
				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id - 1 ) ) = 4.0;
				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id ) ) = -6.0;
				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id + 1 ) ) = 4.0;
				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id + 2 ) ) = -1.0;

				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id - 2 ) + 1 ) = - 1.0;
				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id - 1 ) + 1 ) = 4.0;
				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id ) + 1 ) = -6.0;
				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id + 1 ) + 1 ) = 4.0;
				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id + 2 ) + 1 ) = -1.0;

				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id - 2 ) + 2 ) = - 1.0;
				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id - 1 ) + 2 ) = 4.0;
				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id ) + 2 ) = -6.0;
				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id + 1 ) + 2 ) = 4.0;
				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id + 2 ) + 2 ) = -1.0;


			}



		}


		MATRIX( 3 * ( N_beads - 2 ) , 3 * ( N_beads - 4 ) ) = -1.0;
		MATRIX( 3 * ( N_beads - 2 ) , 3 * ( N_beads - 3 )  ) = 4.0;
		MATRIX( 3 * ( N_beads - 2 ) , 3 * ( N_beads - 2 )  ) = -5.0;
		MATRIX( 3 * ( N_beads - 2 ) , 3 * ( N_beads - 1 )  ) = 2.0;

		MATRIX( 3 * ( N_beads - 2 ) + 1 , 3 * ( N_beads - 4 ) + 1 ) = -1.0;
		MATRIX( 3 * ( N_beads - 2 ) + 1 , 3 * ( N_beads - 3 ) + 1  ) = 4.0;
		MATRIX( 3 * ( N_beads - 2 ) + 1 , 3 * ( N_beads - 2 ) + 1  ) = -5.0;
		MATRIX( 3 * ( N_beads - 2 ) + 1 , 3 * ( N_beads - 1 ) + 1  ) = 2.0;

		MATRIX( 3 * ( N_beads - 2 ) + 2 , 3 * ( N_beads - 4 ) + 2 ) = -1.0;
		MATRIX( 3 * ( N_beads - 2 ) + 2 , 3 * ( N_beads - 3 ) + 2  ) = 4.0;
		MATRIX( 3 * ( N_beads - 2 ) + 2 , 3 * ( N_beads - 2 ) + 2  ) = -5.0;
		MATRIX( 3 * ( N_beads - 2 ) + 2 , 3 * ( N_beads - 1 ) + 2  ) = 2.0;

	}




//bending_matrix

	SparseMatrix<double> b_SPARSE( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList_n_i_mu;
	for( unsigned int row = 0 ; row < 3 * this->numberOfPoints ; row ++ )
	{
		for( unsigned int colum = 0 ; colum < 3 * this->numberOfPoints ; colum ++ )
		{
			if( MATRIX( row , colum ) != 0 )
			{
				tripletList_n_i_mu.push_back(T( row , colum , MATRIX( row , colum ) ));
			}
    		}

    	}

        b_SPARSE.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );
 	this->b_MATRIX = b_SPARSE;


	MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
	for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
	{
		projectionMatrix( index , index ) = 1.0;
	}



	//


	SparseMatrix<double> kronecker( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
        std::vector<T> triplet_kronecker;
	for( unsigned int row = 0 ; row < 3 * this->numberOfPoints ; row ++ )
	{
		for( unsigned int colum = 0 ; colum < 3 * this->numberOfPoints ; colum ++ )
		{
			if( MATRIX( row , colum ) != 0 )
			{
				triplet_kronecker.push_back(T( row , colum , projectionMatrix( row , colum ) ));
			}
    		}

    	}
	this->Kronecker_Delta = kronecker;




/*
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplet_for_sparse;

	for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		Vector3d tangent = this->getTangent2( i );
		tangent = tangent * ( 1.0 / tangent.norm() );
		for( unsigned int j = 0 ; j < 6 ; j ++ )
		{
			int tmp = 3 * i;
			if( j < 3)
			{
				tripletList_n_i_mu.push_back(T( tmp + j , i , - tangent( j ) ));
			}
			else
			{   tripletList_n_i_mu.push_back(T( tmp + j , i , tangent( j - 3 ) ));
			}
		}
	}
	n_i_mu.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );
*/









}



//Bending force calculation using matrixes
void Microtubule::getBendingForces_2( MatrixXd &bendingForce )
{
		//This function enables quick calculation of bending forces using matrix multiplication
		//Curently unused
		if( ( bendingForce.rows() !=  3 * this->numberOfPoints ) || ( ( bendingForce.cols() != 1 ) ) )
		{
			cout<<"( bendingForce.getNN() != 3 * this->numberOfPoints ) || ( ( bendingForce.getNN() != 1 ) ) in Microtubule::getBendingForces( MatDoub bendingForce )"<<endl;
			cout<<this->numberOfPoints <<endl;
			cout<<this->get_dynein_index()<<endl;

			throw("");
		}

		if( this->dydein_index == 20 )
		{
			cout<<"Microtubule ERROR_ID = "<<6466156156135<<endl;
			cout<<" this->dydein_index == 20 "<<endl;
			cout<<"void Microtubule::getBendingForces_2( MatrixXd &bendingForce )"<<endl;
			throw("");
		}
		bendingForce = this->b_MATRIX * this->coordinates * this->kappa / ( this->restDistancePoints * this->restDistancePoints);
}





//Computes bending forces
void Microtubule::getBendingForces( MatrixXd &bendingForce )
{
//https://www.rieger.uni-saarland.de/Paper/supplementary_final_hornak_rieger_2020.pdf

		if( ( bendingForce.rows() !=  3 * this->numberOfPoints ) || ( ( bendingForce.cols() != 1 ) ) )
		{
			cout<<"( bendingForce.getNN() != 3 * this->numberOfPoints ) || ( ( bendingForce.getNN() != 1 ) ) in Microtubule::getBendingForces( MatDoub bendingForce )"<<endl;
			cout<<this->numberOfPoints <<endl;
			cout<<this->get_dynein_index()<<endl;

			throw("");
		}
		if( ( this->getNumberOfPoints() <= 2 )  )
        	{
            		return;
        	}

		if( ( this->getNumberOfPoints() == 3 )  )
        	{

			for( unsigned int index = 0 ; index < this->getNumberOfPoints() ; index ++ )
			{

				Vector3d first = this->getTangent2( 0 );
				Vector3d second = this->getTangent2( 1 );

				double coefficient = second.norm() / ( this->restDistancePoints );
				double treshold = 0.1;
				if( coefficient < treshold )
				{

				}
				else
				{
					coefficient = treshold;
				}
				first = first * ( 1.0 / first.norm() );
				second = second * ( 1.0 / second.norm() );


				coefficient = 1;
				double a_1 = ( - second( 0 ) + first( 0 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;
				double a_2 = ( - second( 1 ) + first( 1 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;
				double a_3 = ( - second( 2 ) + first( 2 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;

				double b_1 = ( first( 0 ) - second( 0 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;
				double b_2 = ( first( 1 ) - second( 1 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;
				double b_3 = ( first( 2 ) - second( 2 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;


				if( index == 0 )										//this is border case of bead of minimum number
				{
					bendingForce( 3 * index +  0 , 0 ) = a_1;
					bendingForce( 3 * index +  1 , 0 ) = a_2;
					bendingForce( 3 * index +  2 , 0 ) = a_3;
				}

				else if( index == 2 )		//this is border case of bead of maximum number
				{
					bendingForce( 3 * index + 0 , 0 ) = b_1;
					bendingForce( 3 * index + 1 , 0 ) = b_2;
					bendingForce( 3 * index + 2 , 0 ) = b_3;

				}

				else if(  index == 1  )
				{
					bendingForce( 3 * index +  0 , 0 ) = ( -1.0 ) * ( a_1 + b_1 );
					bendingForce( 3 * index +  1 , 0 ) = ( -1.0 ) * ( a_2 + b_2 );
					bendingForce( 3 * index +  2 , 0 ) = ( -1.0 ) * ( a_3 + b_3 );

				}
			}
        	}
		else
		{

		    for( unsigned int index = 0 ; index < this->numberOfPoints ; index ++ )
		    {

			if( index == 0 )										//this is border case of bead of minimum number
			{
				Vector3d first( 0.0 , 0.0 , 0.0 );
				Vector3d second( 0.0 , 0.0 , 0.0 );

				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					first( j ) = this->coordinates( 3 * ( index + 1 ) + j , 0 ) - this->coordinates( j , 0 );
					second( j ) =  this->coordinates( 3 * ( index + 2 ) + j , 0 ) - this->coordinates( 3 * ( index + 1 ) + j , 0 );
				}

				double magThisVec = first.norm();
				first = first * ( 1.0 / magThisVec );

				double magNextVec = second.norm();
				second = second * ( 1.0 / magNextVec );

				bendingForce( 3 * index +  0 , 0 ) = ( - second( 0 ) + first( 0 ) * first.dot( second ) ) / magThisVec;
				bendingForce( 3 * index +  1 , 0 ) = ( - second( 1 ) + first( 1 ) * first.dot( second ) ) / magThisVec;
				bendingForce( 3 * index +  2 , 0 ) = ( - second( 2 ) + first( 2 ) * first.dot( second ) ) / magThisVec;
			}

			else if( index == this->numberOfPoints - 1 )		
			{

				Vector3d first( 0.0 , 0.0 , 0.0 );
				Vector3d second( 0.0 , 0.0 , 0.0 );
				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					first( j ) = this->coordinates( 3 * ( index - 1 ) + j , 0 ) - this->coordinates( 3 * ( index - 2 ) + j , 0 );
					second( j ) = this->coordinates( 3 * index + j , 0 ) - this->coordinates( 3 * ( index - 1 ) + j , 0 );
				}

				double magThisVec = first.norm();
				first = first * ( 1.0 / magThisVec );

				double magNextVec = second.norm();
				second = second * ( 1.0 / magNextVec );		

				bendingForce( 3 * index + 0 , 0 ) = ( first( 0 ) - second( 0 ) * first.dot( second ) ) / magNextVec;
				bendingForce( 3 * index + 1 , 0 ) = ( first( 1 ) - second( 1 ) * first.dot( second ) ) / magNextVec;
				bendingForce( 3 * index + 2 , 0 ) = ( first( 2 ) - second( 2 ) * first.dot( second ) ) / magNextVec;

			}
			else
			{
				if( this->numberOfPoints > 3 )
				{

					Vector3d tangent_j = this->getTangent2( index );
                                        Vector3d tangent_j_M_1 = this->getTangent2( index - 1 );

					double tangent_jNorm = tangent_j.norm();
					tangent_j = tangent_j * ( 1.0 / tangent_jNorm );

					double tangent_j_M_1Norm = tangent_j_M_1.norm();
					tangent_j_M_1 = tangent_j_M_1 * ( 1.0 / tangent_j_M_1Norm );


					bendingForce( 3 * index + 0 , 0 ) = ( tangent_j( 0 ) - tangent_j_M_1( 0 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_j_M_1Norm;
					bendingForce( 3 * index + 0 , 0 ) = bendingForce( 3 * index + 0 , 0 ) + ( - tangent_j_M_1( 0 ) + tangent_j( 0 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_jNorm;

					bendingForce( 3 * index + 1 , 0 ) = ( tangent_j( 1 ) - tangent_j_M_1( 1 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_j_M_1Norm;
					bendingForce( 3 * index + 1 , 0 ) = bendingForce( 3 * index + 1 , 0 ) + ( - tangent_j_M_1( 1 ) + tangent_j( 1 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_jNorm;

					bendingForce( 3 * index + 2 , 0 ) = ( tangent_j( 2 ) - tangent_j_M_1( 2 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_j_M_1Norm;
					bendingForce( 3 * index + 2 , 0 ) = bendingForce( 3 * index + 2 , 0 ) + ( - tangent_j_M_1( 2 ) + tangent_j( 2 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_jNorm;


					if ( index == this->numberOfPoints - 2 )
					{

					}

					else
					{
                  				Vector3d first = this->getTangent2( index );
                                                Vector3d second = this->getTangent2( index + 1 );

						double magThisVec = first.norm();
						first = first * ( 1.0 / magThisVec );

						double magNextVec = second.norm();
						second = second * ( 1.0 / magNextVec );		//ATTENTION - here, it will be divided by magNextVec - the last tangent


						bendingForce( 3 * index + 0 , 0 ) = bendingForce( 3 * index + 0 , 0 ) + ( - second( 0 ) + first( 0 ) * tangent_j.dot( second ) ) / magThisVec;
						bendingForce( 3 * index + 1 , 0 ) = bendingForce( 3 * index + 1 , 0 ) + ( - second( 1 ) + first( 1 ) * tangent_j.dot( second ) ) / magThisVec;
						bendingForce( 3 * index + 2 , 0 ) = bendingForce( 3 * index + 2 , 0 ) + ( - second( 2 ) + first( 2 ) * tangent_j.dot( second ) ) / magThisVec;
					}
					if ( index == 1 )
					{

					}
					else
					{
                  				Vector3d first = this->getTangent2( index - 1 );
                                                Vector3d second = this->getTangent2( index - 2 );

						//added to tangent - 2 so, logically, it does not apply for index == 1
						double magThisVec = first.norm();
						first = first * ( 1.0 / magThisVec );

						double magNextVec = second.norm();
						second = second * ( 1.0 / magNextVec );

						bendingForce( 3 * index + 0 , 0 ) = bendingForce( 3 * index + 0 , 0 ) + ( second( 0 ) - first( 0 ) * first.dot( second ) ) / magThisVec;
						bendingForce( 3 * index + 1 , 0 ) = bendingForce( 3 * index + 1 , 0 ) + ( second( 1 ) - first( 1 ) * first.dot( second ) ) / magThisVec;
						bendingForce( 3 * index + 2 , 0 ) = bendingForce( 3 * index + 2 , 0 ) + ( second( 2 ) - first( 2 ) * first.dot( second ) ) / magThisVec;
					}
				}
			}
		    }
		}
		bendingForce = bendingForce * this->kappa;
}








// It returns random forces acting on the whole microtubule
void Microtubule::getRandomForces(  double timeStep , MatrixXd &randomForces )
{
	//The presence of the random forces is controlled by sim_of_Cell::random_force_switch
	std::normal_distribution<> distribution{0,1};
        unsigned int number_of_generator = omp_get_thread_num();

	if( ( randomForces.rows() !=  3 * this->numberOfPoints ) || ( ( randomForces.cols() != 1 ) ) )
	{
		cout<<"( bendingForce.getNN() != 3 * this->numberOfPoints ) || ( ( bendingForce.getNN() != 1 ) ) in Microtubule::getBendingForces( MatDoub bendingForce )"<<endl;
		throw("");
	}

	unsigned int coordinates = 3 * this->numberOfPoints;
	double root = 1.0 * sqrt( 2.0 * sim_of_Cell::boltzmann_constant * sim_of_Cell::Temperature * this->effective_friction / timeStep );

	if( sim_of_Cell::random_force_switch == true )
	{

	}
	else
	{
        	root = 0;
	}
	for( unsigned int i = 0 ; i < coordinates ; i ++ )
	{
		double normal_number = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
		randomForces( i , 0 ) = normal_number * root;
	}
}







//Returns the simplified force of dynein - costant force
Vector3d Microtubule::dynein_force_capture_shrinkage( )
{
	//this function is very simple: just return force acting on the last bead
	Vector3d position_last_bead = this->getPoint( this->getNumberOfPoints() - 1 );
	Vector3d orientace = this->get_IS_position_catching() - position_last_bead;

	orientace = orientace / orientace.norm();
	Vector3d force = orientace * IS_Capture_shrinkage_param::dynein_Force_capture_shrinkage;
	return force;
}











//Mid step algorithm for the propagation of the microtubule
void Microtubule::oneStepMidStepAlgorithm()
{
	//This is one time step to mid-point algorithm
	//Random forces have to be preserved  - same force is used in both halves of step

	//MICROTUBULES
	MatrixXd random_Forces;
	MatrixXd original_Coordinates;
	MatrixXd external_Forces;


	if( this->get_dynein_index() <= 1 )
	{
		original_Coordinates = this->getCoordinates();
		random_Forces = MatrixXd::Zero( 3 * this->getNumberOfPoints() , 1 );
		this->getRandomForces( sim_of_Cell::time_Step , random_Forces );



		//external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
		external_Forces = MatrixXd::Zero( 3 * this->getNumberOfPoints() , 1 );

		//The dynein is added here
		//dynein forces - they apply only when microtubule is cought by dynein motors -
		if( this->get_dynein_index() == 1 )
		{
			//cout<<"this->get_dynein_index() == 1 "<<endl;
			Vector3d dyneinForceVector = this->dynein_force_capture_shrinkage( );

			for( unsigned int i = 0 ; i < 3 ; i ++ )
			{
				external_Forces( 3 * ( this->getNumberOfPoints() - 1 ) + i , 0 ) = dyneinForceVector( i );
			}
		}
		this->oneStepMidStepAlgorithm_1_half( random_Forces , external_Forces );
	}

	//---------------------------------------- SECOND HALF ---------------------------------------------------
	//Allocation of Matrixes of right size

	if( this->get_dynein_index() <= 1 )
	{
		external_Forces = MatrixXd::Zero( 3 * this->getNumberOfPoints() , 1 );
		this->oneStepMidStepAlgorithm_2_half( original_Coordinates , random_Forces , external_Forces );
	}

}




//This is used when calculating movement in the cell
//This is one step of midStep algorithm that includes external forces form MTOC and wall of the cell
void Microtubule::oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )
{

	if( ( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  )
	{
		cout<<"( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" randomForces.rows() = "<<randomForces.rows()<<endl;
		cout<<"randomForces.cols() = "<<randomForces.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		throw("");
	}

	unsigned int rows_tmp = extrenal_Forces.rows();
        unsigned int arg_tmp = this->numberOfPoints * 3;

    if( ( rows_tmp != arg_tmp )  )
	{
                cout<<"ERROR_ID Microtubule6461651685152"<<endl;
		cout<<"( extrenal_Forces.rows() != this->numberOfPoints * 3 ) || ( extrenal_Forces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" extrenal_Forces.rows() = "<<extrenal_Forces.rows()<<endl;
		cout<<"extrenal_Forces.cols() = "<<randomForces.cols()<<endl;
                cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
                cout<<"this->microtubule_id = "<<this->microtubule_id<<endl;
		throw("");
	}

	if( this->numberOfPoints == 0 )
	{
		cout<<"this->numberOfPoints == 0 in Microtubule::oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )"<<endl;
		throw("");
	}
	//Coordinates before step
	MatrixXd original = this->getCoordinates();

	if( this->numberOfPoints == 1 )
	{
		//No projection is needed and also metric forces do not exist.
		MatrixXd V_0 = (  extrenal_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * sim_of_Cell::time_Step_half;
		this->setCoordinates( R_half );
	}
	else
	{
		//Kronecker delta
		MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
		for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
		{
			projectionMatrix( index , index ) = 1.0;
		}

		// Inversin matrix
		MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );
		//Projection matrixes
        	this->getMatrixes( G_uvINV , projectionMatrix );
     		double timeHalfStep =  sim_of_Cell::time_Step_half; /// 2.0

		// ------------------------------------------------------------------------------------------ Inter Microtubule forces
		//Bending Forces
		MatrixXd bendingForce = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
		this->getBendingForces( bendingForce );
		// ------------------------------------------------------------------------------------------ Projections

		MatrixXd forces = bendingForce + extrenal_Forces + randomForces;

		//Projection
		MatrixXd forces_Projection = projectionMatrix * forces;

		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * timeHalfStep;
		this->setCoordinates( R_half );
	}
}


//This is the first half of the mid step algorithm
void Microtubule::oneStepMidStepAlgorithm_1_half_producing_random( MatrixXd& randomForces , MatrixXd& extrenal_Forces )
{
	if( ( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  )
	{
		cout<<"( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" randomForces.rows() = "<<randomForces.rows()<<endl;
		cout<<"randomForces.cols() = "<<randomForces.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		throw("");
	}

	unsigned int rows_tmp = extrenal_Forces.rows();
    	unsigned int arg_tmp = this->numberOfPoints * 3;

    	if( ( rows_tmp != arg_tmp )  )
	{
        	cout<<"ERROR_ID Microtubule6461651685152"<<endl;
		cout<<"( extrenal_Forces.rows() != this->numberOfPoints * 3 ) || ( extrenal_Forces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" extrenal_Forces.rows() = "<<extrenal_Forces.rows()<<endl;
		cout<<"extrenal_Forces.cols() = "<<randomForces.cols()<<endl;
        	cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        	cout<<"this->microtubule_id = "<<this->microtubule_id<<endl;
		throw("");
	}

	if( this->numberOfPoints == 0 )
	{
		cout<<"this->numberOfPoints == 0 in Microtubule::oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )"<<endl;
		throw("");
	}
	//Coordinates before step
	MatrixXd original = this->getCoordinates();

	if( this->numberOfPoints == 1 )
	{
		//No projection is needed and also metric forces do not exist.
		MatrixXd V_0 = (  extrenal_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * sim_of_Cell::time_Step_half;
        	this->setCoordinates( R_half );
	}
	else
	{
		//Kronecker delta
		MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
		for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
		{
			projectionMatrix( index , index ) = 1.0;
		}

		// Inversin matrix
		MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );
        	this->getMatrixes( G_uvINV , projectionMatrix );
          	double timeHalfStep =  sim_of_Cell::time_Step_half; /// 2.0


		// ------------------------------------------------------------------------------------------ Inter Microtubule forces
		//Bending Forces
		MatrixXd bendingForce = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    		this->getBendingForces( bendingForce );


		// ------------------------------------------------------------------------------------------ Projections
		//Projection
		this->getRandomForces(  sim_of_Cell::time_Step , randomForces );

		MatrixXd forces = bendingForce + extrenal_Forces;
		MatrixXd forces_Projection = projectionMatrix * forces  + randomForces;

		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * timeHalfStep;

		this->setCoordinates( R_half );

	}


}

















void Microtubule::oneStepMidStepAlgorithm_2_half( MatrixXd& original_Coordinates , MatrixXd& randomForces , MatrixXd& extrenal_Forces  )
{

	if( this->numberOfPoints == 1 )
	{
		//No projection is needed and also metric forces do not exist.
		MatrixXd V_0 = (  extrenal_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original_Coordinates + V_0 * sim_of_Cell::time_Step;
		this->setCoordinates( R_half );
	}
	else
	{
		//Kronecker delta
		MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
		for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
		{
			projectionMatrix( index , index ) = 1.0;
		}

		// Inversin matrix
		MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );
        	this->getMatrixes( G_uvINV , projectionMatrix );



		// ------------------------------------------------------------------------------------------ Inter Microtubule forces
		//Bending Forces
		MatrixXd bendingForce = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    		this->getBendingForces( bendingForce );

		// ------------------------------------------------------------------------------------------ Projections
		//Projection
		MatrixXd forces = bendingForce + extrenal_Forces;
		MatrixXd forces_Projection = projectionMatrix * forces + randomForces;

		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original_Coordinates + V_0 * sim_of_Cell::time_Step;

        	this->setCoordinates( R_half );

	}
}


//Gets the projection matrix and projects forces
void Microtubule::project_force( MatrixXd& original_force , MatrixXd& projected_forces )
{
	//This function will be used in future simulations
	if( ( original_force.rows() != this->numberOfPoints * 3 ) || ( original_force.cols() != 1 )  )
	{
		cout<<"( original_force.rows() != this->original_force * 3 ) || ( randomForces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" original_force.rows() = "<<original_force.rows()<<endl;
		cout<<"original_force.cols() = "<<original_force.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		cout<<"Microtubule::project_force( MatrixXd& original_force , MatrixXd& projected_forces )"<<endl;
		throw("");
	}
	if( ( projected_forces.rows() != this->numberOfPoints * 3 ) || ( projected_forces.cols() != 1 )  )
	{
		cout<<"( projected_forces.rows() != this->numberOfPoints * 3 ) || ( projected_forces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" projected_forces.rows() = "<<projected_forces.rows()<<endl;
		cout<<"projected_forces.cols() = "<<projected_forces.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		cout<<"Microtubule::project_force( MatrixXd& original_force , MatrixXd& projected_forces )"<<endl;
		throw("");
	}

	MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
	for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
	{
		projectionMatrix( index , index ) = 1.0;
	}

		// Inversin matrix
	MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );
        this->getMatrixes( G_uvINV , projectionMatrix );
	projected_forces = projectionMatrix * original_force;
}





void Microtubule::oneStepMidStepAlgorithm_2_half_projected( MatrixXd& original_Coordinates , MatrixXd& randomForces , MatrixXd& other_Forces  )
{
		if( ( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  )
		{
			cout<<"( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
			cout<<" randomForces.rows() = "<<randomForces.rows()<<endl;
			cout<<"randomForces.cols() = "<<randomForces.cols()<<endl;
			cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
			throw("");
		}
		if( ( other_Forces.rows() != this->numberOfPoints * 3 ) || ( other_Forces.cols() != 1 )  )
		{
			cout<<"( extrenal_Forces.rows() != this->numberOfPoints * 3 ) || ( extrenal_Forces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
			cout<<" extrenal_Forces.rows() = "<<other_Forces.rows()<<endl;
			cout<<"extrenal_Forces.cols() = "<<other_Forces.cols()<<endl;
			cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
			throw("");
		}



		if( this->numberOfPoints == 0 )
		{
			cout<<"this->numberOfPoints == 0 in Microtubule::oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )"<<endl;
			throw("");
		}

	if( this->numberOfPoints == 1 )
	{
		//if this->numberOfPoints == 1 there are no geometrical problems! No projection is needed and also metric forces do not exist.
		//Physically speaking, no changes in inertia caused by configuration

		MatrixXd V_0 = (  other_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original_Coordinates + V_0 * sim_of_Cell::time_Step;
		this->setCoordinates( R_half );
	}
	else
	{
		MatrixXd forces_Projection = randomForces + other_Forces;
		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original_Coordinates + V_0 * sim_of_Cell::time_Step;
        	this->setCoordinates( R_half );

	}
}














//This function propagates the movement of the microtubule with Euler algorithm
void Microtubule::Euler_algorithm( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )
{

    unsigned int rows_tmp = extrenal_Forces.rows();
    unsigned int arg_tmp = this->numberOfPoints * 3;


    if( ( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  )
	{
		cout<<"( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  in Euler_algorithm "<<endl;
		cout<<" randomForces.rows() = "<<randomForces.rows()<<endl;
		cout<<"randomForces.cols() = "<<randomForces.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		throw("");
	}


    if( ( rows_tmp != arg_tmp )  )
	{
        	cout<<"Microtubule ERROR_ID = 13561456464568"<<endl;
		cout<<"( extrenal_Forces.rows() != this->numberOfPoints * 3 ) || ( extrenal_Forces.cols() != 1 )  in Euler_algorithm "<<endl;
		cout<<" extrenal_Forces.rows() = "<<extrenal_Forces.rows()<<endl;
        	cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        	cout<<"this->microtubule_id = "<<this->microtubule_id<<endl;
		throw("");
	}

	if( this->numberOfPoints == 0 )
	{
        	cout<<"ERROR_ID = 132561497813613"<<endl;
		cout<<"this->numberOfPoints == 0 in Microtubule::Euler_algorithm( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )"<<endl;
		throw("");
	}
	//Coordinates before step
	MatrixXd original = this->getCoordinates();

	if( this->numberOfPoints == 1 )
	{
		MatrixXd V_0 = (  extrenal_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * sim_of_Cell::time_Step;
		this->setCoordinates( R_half );
	}
	else
	{
		//Kronecker delta
		MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
		for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
		{
			projectionMatrix( index , index ) = 1.0;
		}

		// Inversin matrix
		MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );
                this->getMatrixes( G_uvINV , projectionMatrix );
		// ------------------------------------------------------------------------------------------ Inter Microtubule forces
		//Bending Forces
		MatrixXd bendingForce = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
		this->getBendingForces( bendingForce );
		MatrixXd forces = bendingForce + extrenal_Forces + randomForces;
		MatrixXd forces_Projection = projectionMatrix * forces; //
		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * sim_of_Cell::time_Step;
	        this->setCoordinates( R_half );
	}
}








void Microtubule::resizeMicrotubule()
{
	//Tangents array contains resized vector
	if( this->getNumberOfPoints() > 2 )
	{
		Vector3d tangents[ this->numberOfPoints - 1 ];
		for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
		{
			Vector3d tangent = this->getTangent2( i );
			tangents[ i ] = tangent * this->restDistancePoints / tangent.norm();
		}
		//Regrow of microtubule
		//Bead 0 stays at the place - in MTOC - other beads are resized
		for( unsigned int i = 1 ; i< this->numberOfPoints ; i ++ )
		{

			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				coordinates( 3 * i + j , 0 ) = coordinates( 3 * ( i - 1 ) + j , 0 ) + tangents[ i - 1 ]( j );
			}
		}
	}

}


//It does the same thing as the previous function - difference is that the first segment has the orientation given by the vector
void Microtubule::resizeMicrotubule( Vector3d orientation )
{
    if( this->getNumberOfPoints() > 2 )
	{
		Vector3d tangents[ this->numberOfPoints - 1 ];

                tangents[ 0 ] = orientation / orientation.norm() * this->restDistancePoints ;
		for( unsigned int i = 1 ; i < this->numberOfPoints - 1 ; i ++ )
		{
			Vector3d tangent = this->getTangent2( i );
			tangents[ i ] = tangent / tangent.norm() * this->restDistancePoints ;
		}

		//Regrow of microtubule
		//Bead 0 stays at the place - in MTOC - other beads are resized
		for( unsigned int i = 1 ; i< this->numberOfPoints ; i ++ )
		{
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				coordinates( 3 * i + j , 0 ) = coordinates( 3 * ( i - 1 ) + j , 0 ) + tangents[ i - 1 ]( j );
			}
		}
	}

}


//Simple resizing of the microtubule
//The first segment has a different size
void Microtubule::resizeMicrotubule_first_segment_different( )
{
	//Tangents array contains resized vector
	if( this->getNumberOfPoints() > 1 )
	{
		Vector3d tangents[ this->numberOfPoints - 1 ];
                Vector3d orientation_1 = this->getTangent2( 0 ) / this->getTangent2( 0 ).norm();


                tangents[ 0 ] = orientation_1 * this->restDistancePoints_first;
		for( unsigned int i = 1 ; i < this->numberOfPoints - 1 ; i ++ )
		{
			Vector3d tangent = this->getTangent2( i );
                        Vector3d orientation = tangent / tangent.norm();
			tangents[ i ] = orientation * this->restDistancePoints;
		}


		//Bead 0 stays at the place - in MTOC - other beads are resized
		for( unsigned int i = 1 ; i< this->numberOfPoints ; i ++ )
		{
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				coordinates( 3 * i + j , 0 ) = coordinates( 3 * ( i - 1 ) + j , 0 ) + tangents[ i - 1 ]( j );
			}
		}
	}
}

//This resized the microtubule with the lengths given by the vector
void Microtubule::resizeMicrotubule_with_different_tangent_lenghts( )
{

	//tohle jen resizne mikrotubuly, neobsahuje odebrani beadu
	//vsechny beady krome prvniho se meni, vcetne posledniho chyceneho v imu syn
	if( this->getNumberOfPoints() > 1 )
	{
		Vector3d tangents[ this->numberOfPoints - 1 ];
        	Vector3d orientation_1 = this->getTangent2( 0 ) / this->getTangent2( 0 ).norm();


        	tangents[ 0 ] = orientation_1 * this->restDistancePoints_first;
		for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
		{
			Vector3d tangent = this->getTangent2( i );
           		Vector3d orientation = tangent / tangent.norm();
			if( std::isnan( orientation.norm() ) )
			{
				cout<<"i = "<<i<<endl;
				cout<<"isnan( orientation.norm() )"<<endl;
				cout<<orientation.norm()<<endl;
				cout<<"tangent.norm() = "<<tangent.norm()<<endl;
				cout<<"  "<<this->getPoint( i )<<endl;
				cout<<"  "<<endl;
				cout<<"  "<<this->getPoint( i + 1 )<<endl;
				cout<<"this->getNumberOfPoints() = "<<this->getNumberOfPoints()<<endl;
				cout<<"void Microtubule::resizeMicrotubule_with_different_tangent_lenghts( )"<<endl;
				throw("");
			}
			if( std::isnan( this->lenght_of_tangents( i , 0 ) ) )
			{
				cout<<"i = "<<i<<endl;
				cout<<"this->lenght_of_tangents( i , 0 ) = "<<this->lenght_of_tangents( i , 0 )<<endl;

				cout<<"PPPPPPPPPPPPPPP"<<endl;
				cout<<"i = "<<i<<endl;
				cout<<"isnan( orientation.norm() )"<<endl;
				cout<<orientation.norm()<<endl;
				cout<<"tangent.norm() = "<<tangent.norm()<<endl;
				cout<<"  "<<this->getPoint( i )<<endl;
				cout<<"  "<<endl;
				cout<<"  "<<this->getPoint( i + 1 )<<endl;
				cout<<"this->getNumberOfPoints() = "<<this->getNumberOfPoints()<<endl;
				cout<<"void Microtubule::resizeMicrotubule_with_different_tangent_lenghts( )"<<endl;
				throw("");
			}



			tangents[ i ] = orientation * this->lenght_of_tangents( i , 0 );
		}
		//Bead 0 stays at the place - in MTOC - other beads are resized
		for( unsigned int i = 1 ; i< this->numberOfPoints ; i ++ )
		{
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				coordinates( 3 * i + j , 0 ) = coordinates( 3 * ( i - 1 ) + j , 0 ) + tangents[ i - 1 ]( j );
			}
		}
	}


}








void Microtubule::Substract_one_Bead()
{
	if( this->numberOfPoints == 1  )
	{
		// this means that there is important misuse in calling of the function
		cout<<"this->numberOfPoints == 1  and I call Microtubule::Substract_one_Bead()"<<endl;
		throw("");
	}

    if( this->numberOfPoints > 2 )
    {
        this->numberOfPoints = this->numberOfPoints - 1;
        if( this->numberOfPoints == 1 )
        {
            cout<<".............SUBSTRACTION....."<<endl;
        }

        MatrixXd tmp = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
        for( unsigned int i = 0 ; i < 3 * this->numberOfPoints ; i ++ )
        {
            tmp( i , 0 ) = this->coordinates( i , 0 );
        }

        this->coordinates = tmp;
    }
    else
    {
        this->set_dynein_index( 2 );
    }
}


//Simply removes the last bead of the microtubule
void Microtubule::Substract_one_Bead_simple(  )
{
    if( this->numberOfPoints > 2 )
    {
        MatrixXd tmp = MatrixXd::Zero( 3 * ( this->numberOfPoints - 1 ) , 1 );
        //Adjust coordinates
        for( unsigned int i = 0 ; i < 3 * ( this->numberOfPoints - 1 ) ; i ++ )
        {
            tmp( i , 0 ) = this->coordinates( i , 0 );
        }
        //Adjust the number of beads
        this->setNumberOfPoints( this->getNumberOfPoints() - 1 );
	this->setCoordinates( tmp );
    }
    else
    {
        return;
    }
}





//Places the last bead of the microtubule to the position, where microtubule depolymerizes
void Microtubule::add_one_final_Bead_catching_position()
{
        this->numberOfPoints = this->numberOfPoints + 1;
        if( this->numberOfPoints == 1 )
        {
            cout<<".............SUBSTRACTION....."<<endl;
        }
        //copiing old coordinates
        MatrixXd tmp = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
        //coordinates are one long vector: iteration through just one dimension is possible
        for( unsigned int i = 0 ; i < 3 * ( this->numberOfPoints - 1 ) ; i ++ )
        {
            tmp( i , 0 ) = this->coordinates( i , 0 );
        }

        for( unsigned int i = 0 ; i < 3 ; i ++ )
        {
            tmp( 3 * ( this->numberOfPoints - 1 ) + i , 0 ) = this->IS_position_catching( i );
        }

        this->coordinates = tmp;

}














//This prints the beads into terminal - basic control
void Microtubule::print_Microtubule()
{
	//Creation of name of microtubule using sprintf

	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		cout<<" i = "<<i<<endl;
		cout<<this->getPoint( i )<<endl;
	}

}

 //This prints the tangents of microtubules into terminal - basic control
void Microtubule::print_Microtubule_Tangent()
{
    for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		cout<<" i = "<<i<<endl;
		cout<<this->getTangent2( i ).norm()<<endl;
	}
	cout<<"resting_distance = "<<this->getRestDist()<<endl;

}



//It determines whether the point is in the cell 
bool Microtubule::confirm_inner_ellipsoid( Vector3d position , double a_axis , double b_axis )
{
	//It is used in the constructors to heck the position of the beads
	if( position( 2 ) > 0 )
	{
		double confirm_ellipsoid_value = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) / ( a_axis * a_axis );
		confirm_ellipsoid_value = confirm_ellipsoid_value + ( position( 2 ) * position( 2 ) ) / ( b_axis * b_axis );
		if( confirm_ellipsoid_value < 1.0 )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	if( position( 2 ) <= 0 )
	{
		//Long MTs touch the membrane just with the plus end
		b_axis = b_axis - 1.0e-6; 
		double confirm_ellipsoid_value = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) / ( a_axis * a_axis );
		confirm_ellipsoid_value = confirm_ellipsoid_value + ( position( 2 ) * position( 2 ) ) / ( b_axis * b_axis );
		if( confirm_ellipsoid_value < 1.0 )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return false;

}


//The force acting on a point on the segment is transmitted to two beads at the end of segments
void Microtubule::distribute_force_on_beads( MatrixXd& force_on_microtubule , Vector3d force , unsigned int bead_segment , double t )
{
    //No sense to compute this function with zero force in this case anyway
    if( force.norm() < 1e-16 )
    {
        return;
    }

    if( bead_segment >= this->numberOfPoints )
    {
        cout<<"bead_segment >= this->numberOfPoints"<<endl;
        cout<<"Microtubule ERROR_ID = 646161654346468"<<endl;
        cout<<"bead_segment = "<<bead_segment<<endl;
        cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
        cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        throw("");
    }
        Vector3d tangent = this->getTangent2( bead_segment );
        //calculation of the cosinus
        double cosinus = tangent.dot( force ) / ( tangent.norm() * force.norm() );

        //PARALLEL FORCES
        Vector3d parallel_forces = ( force.norm() * cosinus ) * ( tangent / tangent.norm() );


        for( unsigned int i = 0 ; i < 3 ; i ++ )
        {
            force_on_microtubule( 3 * bead_segment + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + parallel_forces( i ) / 2.0;
            force_on_microtubule( 3 * ( bead_segment + 1 ) + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + parallel_forces( i ) / 2.0;
        }

        //PERPENDICULAR FORCES
        Vector3d perpendicullar_forces = force - parallel_forces;
        for( unsigned int i = 0 ; i < 3 ; i ++ )
        {
            force_on_microtubule( 3 * bead_segment + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + perpendicullar_forces( i ) / 2.0;
            force_on_microtubule( 3 * ( bead_segment + 1 ) + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + perpendicullar_forces( i ) / 2.0;
        }

        //MOMENT FORCES
        double ratio = abs( 0.5 - t );
        Vector3d moment_forces = perpendicullar_forces * ratio;

        if( t < 0.5 )
        {
            for( unsigned int i = 0 ; i < 3 ; i ++ )
            {
                force_on_microtubule( 3 * bead_segment + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + moment_forces( i ) / 2.0;
                force_on_microtubule( 3 * ( bead_segment + 1 ) + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) - moment_forces( i ) / 2.0;
            }
        }
        else
        {
            for( unsigned int i = 0 ; i < 3 ; i ++ )
            {
                force_on_microtubule( 3 * bead_segment + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) - moment_forces( i ) / 2.0;
                force_on_microtubule( 3 * ( bead_segment + 1 ) + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + moment_forces( i ) / 2.0;
            }
        }

}






//This adds one pair - anchor point and abscissa to microtubule
void Microtubule::add_pair( std::pair < Vector3d ,double  > point_abscissa )
{
    this->Dynein_motors_2.push_back( point_abscissa );
}


//Returns all dyneins acting on microtubules and erase them
std::vector< Vector3d  > Microtubule::get_Dynein_points_and_erase()
{
    if( this->get_dynein_index() != 9 )
    {
        cout<<"this->get_dynein_index() != 9"<<endl;
        cout<<"std::vector< Vector3d  > Microtubule::get_Dynein_points_and_erase()"<<endl;
        cout<<"ERROR_ID = 6835461681364813"<<endl;
        throw("");
    }

    std::vector< Vector3d  > point_coordinates;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d anchor_position = std::get<0>( pair_tmp );
        //Vector3d anchor_position = pair_tmp[ 0 ];
        point_coordinates.push_back( anchor_position );
    }
    this->Dynein_motors_2.clear();
    return point_coordinates;
}


//It does the same a the previous function
std::vector<Vector3d> Microtubule::get_dynein_points_in_IS_and_erase()
{
    if( ( this->get_dynein_index() != 20 ) && ( this->get_dynein_index() != 40 )  )
    {
        cout<<"this->get_dynein_index() != 20"<<endl;
        cout<<"std::vector< Vector3d  > Microtubule::get_dynein_points_in_IS_and_erase()"<<endl;
        unsigned int ERROR_ID = 781646134;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }

    std::vector< Vector3d  > point_coordinates;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d anchor_position = std::get<0>( pair_tmp );
        point_coordinates.push_back( anchor_position );
    }
    this->Dynein_motors_2.clear();
    return point_coordinates;

}


//Returns all the dyneins acting on microtubule
std::vector< Vector3d  > Microtubule::get_Dynein_points_without_erasing()
{

    if( this->get_dynein_index() < 3  )
    {
        cout<<"this->get_dynein_index() < 3"<<endl;
        cout<<"std::vector< Vector3d  > Microtubule::get_Dynein_points_without_erasing()"<<endl;
        cout<<"ERROR_ID = 97964651643515"<<endl;
        throw("");
    }

    std::vector< Vector3d  > point_coordinates;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d anchor_position = std::get<0>( pair_tmp );
        point_coordinates.push_back( anchor_position );
    }
    return point_coordinates;
}






//Returns force of dynein acting on microtubule
Vector3d Microtubule::force_real_dynein_one_pair(  std::pair < Vector3d , double > pair_tmp   )
{
    Vector3d force( 0.0 , 0.0 , 0.0 );
    Vector3d anchor_position = std::get<0>( pair_tmp );
    double abscissa = std::get<1>( pair_tmp );
    Vector3d point_of_attachment = get_attachment_point_according_to_abscissa( abscissa );
    Vector3d ancher_atachment = anchor_position - point_of_attachment;

    double distance_ancher_atachment = ancher_atachment.norm();

    if( distance_ancher_atachment < Dynein_real::L_0 )
    {

    }
    else
    {
	//cout<<"distance_ancher_atachment = "<<distance_ancher_atachment<<endl;
        double distance_spring = distance_ancher_atachment - Dynein_real::L_0;
        double force_absolute_value = Dynein_real::alpha * distance_spring;
        force = force_absolute_value * ancher_atachment / ancher_atachment.norm();

    }
    return force;


}






//This function performs the steppping and detachment of dynein on the microtubule
std::vector< Vector3d > Microtubule::stepping_detach_real_dynein_abscissa_projection( )
{
    //This function does the stepping of dynein
    //https://www.rieger.uni-saarland.de/Paper/supplementary_final_hornak_rieger_2020.pdf

    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();


    std::vector< Vector3d > return_vectors;
    std::vector< std::pair < Vector3d , double > > new_vectors;


    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d ,double > pair_tmp = this->Dynein_motors_2[ pair_index ];
        Vector3d anchor_position = std::get<0>( pair_tmp );
        double abscissa = std::get<1>( pair_tmp );

        Vector3d force_one_pair = this->force_real_dynein_one_pair(  pair_tmp  );
	double force = force_one_pair.norm();

        unsigned int lower_bead_index = 0;
        if( abscissa < this->getTangent2( 0 ).norm() )
        {
            lower_bead_index = 0;
        }
        else
        {
            unsigned int lower_bead_index;
            try
            {
                lower_bead_index = this->get_index_according_to_abscissa( abscissa );
            }
            catch( unsigned int error_id )
            {
                cout<<"std::vector< Vector3d > Microtubule::stepping_detach_real_dynein_abscissa_projection( )"<<endl;
                throw("");
            }
        }


        Vector3d tangent( 0.0 , 0.0 , 0.0 );
        try
        {
            if( abscissa < this->get_lenght_of_microtubule() )
            {
                tangent = this->getTangent2( lower_bead_index );
            }
            else if( abscissa == this->get_lenght_of_microtubule() )
            {
                tangent = ( -1.0 ) * this->getTangent2( this->getNumberOfPoints() - 2 );
            }
        }
        catch ( string e)
        {
            cout<<"error in std::vector< Vector3d > Microtubule::stepping_detach_real_dynein_abscissa_projection( )"<<endl;
        }

        double detach_prob = Dynein_real::detachment_probability_per_step( force );//NARAZNIK

	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        if( probability < detach_prob )
        {
            //Add the vector to the dynein on surface
            return_vectors.push_back( anchor_position );
            continue;
        }

        if( force == 0 )
        {
 	    double probability_2 = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
            if( probability_2 < Dynein_real::forward_stepping_probability_per_step )
            {
                abscissa = abscissa - Dynein_real::step;
            }
        }
        else
        {
            double cosinus = tangent.dot( force_one_pair ) / ( force_one_pair.norm() * tangent.norm() );
            if( cosinus < 0 )
            {
                //It means that the force on the atachment forces the atachment point to go in the supposed direction ( - )
                //It will add the force to upper bead
                double probability_2 = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                if( probability_2 < Dynein_real::forward_stepping_probability_per_step )
                {
                    abscissa = abscissa - Dynein_real::step;
                }
            }
            else
            {
                //The force aim to the + end            
                if( force < Dynein_real::F_stall )
                {
                    double prav_tmp = Dynein_real::backward_stepp_prob_smaller_than_Stall_f_per_step( force );
                    //double probability_2 = rand_x( 0.0 ,  1.0 );
                    double probability_2 = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                    if( probability_2 < prav_tmp )
                    {
                        abscissa = abscissa - Dynein_real::step;
                    }

                }
                else
                {
                    //It goes to the + end                 
                    double tmp = Dynein_real::backward_prob_force_bigger_then_Stall_force_per_step;
                    //double probability_2 = rand_x( 0.0 ,  1.0 );
                    double probability_2 = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                    if( probability_2 < tmp )
                    {
                        abscissa = abscissa + Dynein_real::step;
                    }

                }
            }
        }
        //The first segment is a part of the MTOC
        if( abscissa <= this->getTangent2( 0 ).norm()  )
        {
            return_vectors.push_back( anchor_position );
        }
        else if( abscissa > this->get_lenght_of_microtubule() )
        {
            //The substraction is because of numerical impreccision when simulating the polymers
            //It can reach the tangent that does not exist
            return_vectors.push_back( anchor_position );
        }
        else
        {
            std::pair < Vector3d ,double > pair_2( anchor_position , abscissa );
            new_vectors.push_back( pair_2 );
        }

    }

    this->Dynein_motors_2 = new_vectors;
    if( new_vectors.size() == 0 )
    {
        if( this->get_dynein_index() == 9 )
        {
            this->set_dynein_index( 0 );
        }
    }

    return return_vectors;
}





//Check whether dyneins are not detached due to the depolymerization
std::vector< Vector3d > Microtubule::detachmet_due_to_shrinkage_of_MT_9( )
{

    if( this->get_dynein_index() != 9 )
    {
	cout<<"this->get_dynein_index() != 9" <<endl;
        cout<<"std::vector< Vector3d > Microtubule::detachmet_due_to_shrinkage_of_MT_9( )"<<endl;
        cout<<"ERROR_ID = 133243251564536"<<endl;
	throw("");
    }		


    std::vector< Vector3d > return_vectors;
    std::vector< std::pair < Vector3d , double > > new_vectors;

    //if the microtubule is shorter than the abscissa, dynein detaches 
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d ,double > pair_tmp = this->Dynein_motors_2[ pair_index ];
        Vector3d anchor_position = std::get< 0 >( pair_tmp );
        double abscissa = std::get< 1 >( pair_tmp );

	if( abscissa < this->getTangent2( 0 ).norm() )
	{
		return_vectors.push_back( anchor_position );
	}
	else if( abscissa >= this->get_lenght_of_microtubule() )
	{
		return_vectors.push_back( anchor_position );
	}
	else
	{
		new_vectors.push_back( pair_tmp );
	}

    }

    this->Dynein_motors_2 = new_vectors;
    if( new_vectors.size() == 0 )
    {
        if( this->get_dynein_index() == 9 )
        {
            this->set_dynein_index( 0 );
        }
    }

    return return_vectors;
}















//Matrix contain all the forces from the dynein
void Microtubule::force_dynein_abscissa_real_dynein( MatrixXd& force_dynein )
{
    //In this function I will just iterate throught all dynein and calculate all the forces

    if( ( force_dynein.rows() != 3 * this->getNumberOfPoints() ) || ( force_dynein.cols() != 1 )  )
    {
        cout<<"( force_dynein.rows() != 3 * this->getNumberOfPoints() ) || ( force_dynein.cols() != 1 )"<<endl;
        cout<<"Microtubule::force_dynein_abscissa_real_dynein( MatrixXd& force_dynein )"<<endl;
        cout<<"Microtubule ERROR_ID = 3516145415461115626"<<endl;
        throw("");
    }


    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d force = this->force_real_dynein_one_pair( pair_tmp );

        Vector3d anchor_position = std::get<0>( pair_tmp );
        double abscissa = std::get<1>( pair_tmp );

        unsigned int index_lower_bead = this->get_index_according_to_abscissa( abscissa );


        if( index_lower_bead < this->getNumberOfPoints() - 1 )
        {

            Vector3d tangent = this->getTangent2( index_lower_bead );
            double abscissa_minus_lower_bead = abscissa - get_distance_to_lower_bead_with_index( index_lower_bead );
            double t_parameter = abscissa_minus_lower_bead / tangent.norm();
            this->distribute_force_on_beads(  force_dynein , force , index_lower_bead , t_parameter );
        }
        else if( index_lower_bead == this->getNumberOfPoints() - 1 )
        {
            double t_parameter = 1.0;
            this->distribute_force_on_beads(  force_dynein , force ,  index_lower_bead - 1 , t_parameter );
        }
    }
}






//Returns the number of segment(or lower bead), in which abscissa is located
unsigned int Microtubule::get_index_according_to_abscissa( double abscissa )
{

    if( abscissa > this->get_lenght_of_microtubule() )
    {
        cout<<"............................................"<<endl;
        cout<<this->get_polygon_number()<<endl;
        cout<<this->getID()<<endl;
        cout<<"dynein index = "<<this->get_dynein_index()<<endl;
        cout<<"abscissa >= this->get_lenght_of_microtubule()"<<endl;
        double diffference = abscissa - this->get_lenght_of_microtubule();
        printf( "%.50f", diffference );
        cout<<endl;


        cout<<"abscissa = "<<endl;
        printf( "%.50f\n", abscissa );
        cout<<"lenght = "<<endl;
        printf( "%.50f\n", this->get_lenght_of_microtubule() );
        cout<<"time = "<<Cell_parametres::time<<endl;


        unsigned int ERROR_ID = 81346;
        cout<<"unsigned int Microtubule::get_index_according_to_abscissa( double abscissa )"<<endl;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }
    else if( abscissa == this->get_lenght_of_microtubule() )
    {
        return ( this->numberOfPoints - 1);
    }
    else
    {



        double lower_length = 0;
        double upper_length = 0;
        for( unsigned int segment_id = 0 ; segment_id < this->getNumberOfPoints() - 1 ; segment_id ++ )
        {
            double length_of_segment = this->getTangent2( segment_id ).norm();
            //cout<<"length_of_segment = "<<length_of_segment<<endl;
            upper_length = upper_length + length_of_segment;
            if( ( abscissa < upper_length ) && ( abscissa >= lower_length ) )
            {
                return segment_id;
            }

            lower_length = lower_length + length_of_segment;
        }
    }
    return 0;
}


//Returns attachment point according to the abscissa
Vector3d Microtubule::get_attachment_point_according_to_abscissa( double abscissa )
{
    if( abscissa == this->get_lenght_of_microtubule() )
    {
        return this->getPoint( this->numberOfPoints - 1 );
    }
    else
    {

        //cout<<"Vector3d Microtubule::get_attachment_point_according_to_abscissa( double abscissa ) 1 "<<endl;
        unsigned int lower_index;
        try
        {
            lower_index = this->get_index_according_to_abscissa( abscissa );
        }
        catch( unsigned int error_id )
        {
            cout<<"std::vector< Vector3d > Microtubule::get_attachment_point_according_to_abscissa( )"<<endl;
            throw("");
        }

        double lower_length = 0;
        for( unsigned int segment = 0 ; segment < lower_index ; segment ++ )
        {
            lower_length = lower_length + this->getTangent2( segment ).norm();
        }

        double abscissa_minus_lower_bead = abscissa - lower_length;
        Vector3d micro_bead = this->getPoint( lower_index );
        Vector3d tangent = this->getTangent2( lower_index );
        Vector3d attachment_point = micro_bead + tangent / tangent.norm() * abscissa_minus_lower_bead;
        return attachment_point;
    }
}



//The lenght of all tangent are saved 
void Microtubule::set_lenght_of_tangents( string argument )
{
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int counter = 0 ; counter < this->numberOfPoints - 1 ; counter ++ )
    {
        Vector3d tangent;
        try
        {
            tangent = this->getTangent2( counter );
        }
        catch( int e )
        {
            cout<<"exception"<<endl;
	    cout<<"void Microtubule::set_lenght_of_tangents( string argument )"<<endl;
	    cout<<"this->cell_id = "<<this->cell_id<<endl;
            throw("");
        }
	if( this->getTangent2( counter ).norm() < 1e-9 )
	{
		cout<<"this->getTangent2( counter ).norm() < 1e-9"<<endl;
		cout<<this->getTangent2( counter ).norm()<<endl;
		cout<<"void Microtubule::set_lenght_of_tangents( string argument )"<<endl;
		cout<<"this->cell_id = "<<this->cell_id<<endl;
		cout<<"argument = "<<argument<<endl;
		cout<<"growth_index = "<<this->growth_index<<endl;
		cout<<"dydein_index = "<<this->dydein_index<<endl;
		cout<<"this->getNumberOfPoints() = "<<this->getNumberOfPoints()<<endl;
		cout<<this->coordinates<<endl;
		throw("");
	}
        this->lenght_of_tangents( counter , 0 ) = this->getTangent2( counter ).norm();
    }
}


//Control of the microtubule length and the stepping of dynein
std::vector<Vector3d> Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule()
{
    //currently unused
    //first, I have to resize micros with old lenghts of segments
    //second, motors will perform stepping
    //third, micro is resized with depolimerization
    //fourth, control if motors willl not be kicked by depolimerization: abscissa longer than lenght of micro

    if( this->get_dynein_index() != 20 )
    {
        cout<<"Microtubule::control_and_resizing_of_depolimerizing_microtubule()"<<endl;
        cout<<" this->get_dynein_index() != 20"<<endl;
        unsigned int ERROR_ID = 5614681;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }




    //First
    this->resizeMicrotubule_with_different_tangent_lenghts();

    //Second
    std::vector<Vector3d> detached_motors = this->stepping_detach_real_dynein_abscissa_projection();




    if( this->getNumberOfPoints() > 3 )
    {
        if( this->getTangent2( this->numberOfPoints - 2 ).norm() < this->restDistancePoints / 2.0   )
        {
            this->numberOfPoints = this->numberOfPoints - 1;
            MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
            for( unsigned int counter = 0 ; counter < 3 * ( this->numberOfPoints - 1 ) ; counter ++ )
            {
                tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
            }

            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                tmp_new( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
            }
            this->coordinates = tmp_new;


        }
        else
        {
            //setting the last point
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                this->coordinates( 3 * ( this->getNumberOfPoints() - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );

            }
        }
        this->set_lenght_of_tangents("std::vector<Vector3d> Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule()");
    }


    else if( this->getNumberOfPoints() == 3 )
    {
            //setting the last point
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            this->coordinates( 3 * ( this->getNumberOfPoints() - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
        }

        this->set_lenght_of_tangents("std::vector<Vector3d> Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule()");





    }


    std::vector<Vector3d> vectors_to_be_returned_to_IS;

    for( unsigned int counter = 0 ; counter < detached_motors.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( detached_motors[ counter ] );
    }

    std::vector<Vector3d> kicked_due_to_repolarization = this->control_motor_detachment_with_depolimerization();
    for( unsigned int counter = 0 ; counter < kicked_due_to_repolarization.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( kicked_due_to_repolarization[ counter ] );
    }


    if( this->Dynein_motors_2.size() == 0 )
    {
 
    }
    if( this->getNumberOfPoints() <= 2 )
    {
        this->set_dynein_index( 0 );
    }


    return vectors_to_be_returned_to_IS;


}




std::vector<Vector3d> Microtubule::get_dynein_points_Cortical_Sliding()
{
    if( this->get_dynein_index() != 9 )
    {
        cout<<"this->get_dynein_index() != 9"<<endl;
        cout<<"std::vector< Vector3d  > Microtubule::get_dynein_points_Cortical_Sliding()"<<endl;
        unsigned int ERROR_ID = 7816134;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }

    std::vector< Vector3d  > point_coordinates;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d anchor_position = std::get<0>( pair_tmp );
        point_coordinates.push_back( anchor_position );
    }
    return point_coordinates;

}


//Checks the last segment of the microtubule and resizes it
std::vector<Vector3d>  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2()
{

    if( ( this->get_dynein_index() != 20 ) && ( this->get_dynein_index() != 40 ) )
    {
        cout<<"Microtubule::control_and_resizing_of_depolimerizing_microtubule()"<<endl;
        cout<<" this->get_dynein_index() != 20"<<endl;
        unsigned int ERROR_ID = 5614681;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }

    this->resizeMicrotubule_with_different_tangent_lenghts();
    std::vector<Vector3d> detached_motors = this->stepping_detach_real_dynein_abscissa_projection();


    if( this->numberOfPoints > 3 )
    {
       //check, whether it depolimerized
       if( this->getTangent2( this->numberOfPoints - 2 ).norm() < 0.9 * sim_of_Cell::resting_distance  )
       {

	    double lenght_micro_outside_MTOC_2 = this->get_lenght_of_microtubule_outside_MTOC();
	    unsigned int number_of_points_outside_MTOC_2 = this->getNumberOfPoints() - 2;
	    double a_a = abs( lenght_micro_outside_MTOC_2 / ( double ) number_of_points_outside_MTOC_2 - sim_of_Cell::resting_distance  );
	    double b_b = abs( lenght_micro_outside_MTOC_2 / ( double ) ( number_of_points_outside_MTOC_2 - 1.0 ) - sim_of_Cell::resting_distance );

	   if( b_b <= a_a )
	   {
		std::vector< Vector3d > orientations;
	    	for( unsigned int point_id = 1 ; point_id < this->getNumberOfPoints() - 1 ; point_id ++ )
	    	{
			orientations.push_back( this->getTangent2( point_id ) / this->getTangent2( point_id ).norm() );
	    	}
	    	double lenght_micro_outside_MTOC = this->get_lenght_of_microtubule_outside_MTOC();

	    	//One bead is removed
	   	this->numberOfPoints = this->numberOfPoints - 1;
	    	unsigned int number_of_points_outside_MTOC = this->getNumberOfPoints() - 2;
		if( this->getNumberOfPoints() - 3 > 0 )
		{
	    		this->restDistancePoints = lenght_micro_outside_MTOC / ( double ) ( number_of_points_outside_MTOC );
		}
		else
		{
			this->restDistancePoints = lenght_micro_outside_MTOC;
	    	}

            	MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    	    	for( unsigned int counter = 0 ; counter < 6 ; counter ++ )
    	    	{
			tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
    	    	}
    	   	for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
    	   	{
			Vector3d orient = this->restDistancePoints * orientations[ i - 2 ];
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				tmp_new( 3 * i + j , 0 ) = tmp_new( 3 * ( i - 1 ) + j , 0 ) + orient( j );
			}
    	   	}


    	  	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	  	{
			tmp_new( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	  	}

    	  	this->coordinates = tmp_new;
	   }
           else
	   {
		std::vector< Vector3d > orientations;
	    	for( unsigned int point_id = 1 ; point_id < this->getNumberOfPoints() - 1 ; point_id ++ )
	    	{
			orientations.push_back( this->getTangent2( point_id ) / this->getTangent2( point_id ).norm() );
	    	}
	    	double lenght_micro_outside_MTOC = this->get_lenght_of_microtubule_outside_MTOC();

	    	//The number of beads stays the same
	   	this->numberOfPoints = this->numberOfPoints;
	    	unsigned int number_of_points_outside_MTOC = this->getNumberOfPoints() - 2;
	    	//It is already changed for shorter microtubule
	    	this->restDistancePoints = lenght_micro_outside_MTOC / ( double ) ( number_of_points_outside_MTOC );

            	MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    	    	for( unsigned int counter = 0 ; counter < 6 ; counter ++ )
    	    	{
			tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
    	    	}

    	   	for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
    	   	{
			Vector3d orient = this->restDistancePoints * orientations[ i - 2 ];
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				tmp_new( 3 * i + j , 0 ) = tmp_new( 3 * ( i - 1 ) + j , 0 ) + orient( j );
			}
    	   	}

    	  	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	  	{
			tmp_new( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	  	}
    	  	this->coordinates = tmp_new;


    	       for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	       {
		    this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	       }
	   }


       }
       else
       {
    	   for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	   {
		this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	   }

       }
   }
   else
   {

       for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
       {
	   this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
       }

   }



    this->set_lenght_of_tangents("std::vector<Vector3d>  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2()");




    std::vector<Vector3d> vectors_to_be_returned_to_IS;

    for( unsigned int counter = 0 ; counter < detached_motors.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( detached_motors[ counter ] );
    }

    std::vector<Vector3d> kicked_due_to_repolarization = this->control_motor_detachment_with_depolimerization();


    for( unsigned int counter = 0 ; counter < kicked_due_to_repolarization.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( kicked_due_to_repolarization[ counter ] );
    }


    if( this->Dynein_motors_2.size() == 0 )
    {
        //this->set_dynein_index( 0 );
    }
    if( this->getNumberOfPoints() <= 2 )
    {
        this->set_dynein_index( 0 );
    }

    return vectors_to_be_returned_to_IS;



}










//Controls the last segment, resizes microtubule and checks the detachment of dyneins
std::vector<Vector3d>  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple()
{

    if( ( this->get_dynein_index() != 20 ) && ( this->get_dynein_index() != 40 ) )
    {
        cout<<"Microtubule::control_and_resizing_of_depolimerizing_microtubule()"<<endl;
        cout<<" this->get_dynein_index() != 20"<<endl;
        unsigned int ERROR_ID = 5614681;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }
    //Resizing
    this->resizeMicrotubule_with_different_tangent_lenghts();
    //Stepping
    std::vector<Vector3d> detached_motors = this->stepping_detach_real_dynein_abscissa_projection();
    if( this->numberOfPoints > 3 )
    {

       if( this->getTangent2( this->numberOfPoints - 2 ).norm() < 0.4 * sim_of_Cell::resting_distance  )
       {

	   	this->numberOfPoints = this->numberOfPoints - 1;
            	MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    	    	for( unsigned int counter = 0 ; counter < 3 * this->numberOfPoints ; counter ++ )
    	    	{
			tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
    	    	}

    	  	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	  	{
			tmp_new( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	  	}

    	  	this->coordinates = tmp_new;

       }
       else
       {
    	   for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	   {
		this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	   }

       }
   }
   else
   {

       for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
       {
	   this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
       }

       if( this->get_last_Tangent( ).norm() < IS_Capture_shrinkage_param::depolimerization_treshold )
       {
		std::vector<Vector3d> vectors_to_be_returned_to_IS  = this->detach_from_IS_2();
		this->set_growth_index( 0 );
    		this->set_lenght_of_tangents("std::vector<Vector3d>  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple()");
		this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );
		for( unsigned int index_id = 0 ; index_id < detached_motors.size() ; index_id ++ )
		{
			vectors_to_be_returned_to_IS.push_back( detached_motors[ index_id ] );
		}
    		return vectors_to_be_returned_to_IS;
       } 		
   }

    this->set_lenght_of_tangents("std::vector<Vector3d>  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple()");


    std::vector<Vector3d> vectors_to_be_returned_to_IS;

    for( unsigned int counter = 0 ; counter < detached_motors.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( detached_motors[ counter ] );
    }

    std::vector<Vector3d> kicked_due_to_repolarization = this->control_motor_detachment_with_depolimerization();
    for( unsigned int counter = 0 ; counter < kicked_due_to_repolarization.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( kicked_due_to_repolarization[ counter ] );
    }


    if( this->Dynein_motors_2.size() == 0 )
    {
        //this->set_dynein_index( 0 );
    }

    return vectors_to_be_returned_to_IS;

}





//Resizing of the microtubule for the microtubule attached to capture-shrinkage dynein
void  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple_after_catching()
{

    if( ( this->get_dynein_index() != 20 ) && ( this->get_dynein_index() != 40 ) )
    {
        cout<<"Microtubule::control_and_resizing_of_depolimerizing_microtubule()"<<endl;
        cout<<" this->get_dynein_index() != 20"<<endl;
        unsigned int ERROR_ID = 5614681;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }

    this->resizeMicrotubule_with_different_tangent_lenghts();
    std::vector<Vector3d> detached_motors = this->stepping_detach_real_dynein_abscissa_projection();


    if( this->numberOfPoints > 3 )
    {

       if( this->getTangent2( this->numberOfPoints - 2 ).norm() < 0.4 * sim_of_Cell::resting_distance  )
       {

	   	this->numberOfPoints = this->numberOfPoints - 1;
            	MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    	    	for( unsigned int counter = 0 ; counter < 3 * this->numberOfPoints ; counter ++ )
    	    	{
			tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
    	    	}

    	  	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	  	{
			tmp_new( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	  	}

    	  	this->coordinates = tmp_new;

       }
       else
       {
    	   for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	   {
		this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	   }

       }
   }
   else
   {

       for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
       {
	   this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
       }

   }



    this->set_lenght_of_tangents("void  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_simple_after_catching()");
}























//Controls the the depolymerizing microtubule and resize it when necessary 
void  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0()
{

    if( ( this->get_dynein_index() != 0 ) )
    {
        cout<<"Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0()"<<endl;
        cout<<" this->get_dynein_index() != 0"<<endl;
        unsigned int ERROR_ID = 456656113;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }

    if( ( this->getNumberOfPoints() <= 2 ) )
    {
        cout<<"Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0()"<<endl;
        cout<<"this->getNumberOfPoints() <= 2"<<endl;
        unsigned int ERROR_ID = 13212323;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }

    double lenght_1 = this->get_lenght_of_microtubule_outside_MTOC();


    this->resizeMicrotubule_with_different_tangent_lenghts(  );

   unsigned int counter = 0;


    if( this->numberOfPoints > 3 )
    {
	//One bead can be removed
       if( this->getTangent2( this->numberOfPoints - 2 ).norm() < 0.5 * sim_of_Cell::resting_distance  ) //0.5
       {
	    counter = counter + 1;
	    double lenght_micro_outside_MTOC_2 = this->get_lenght_of_microtubule_outside_MTOC();
	    unsigned int number_of_points_outside_MTOC_2 = this->getNumberOfPoints() - 2;
	    double a_a = abs( lenght_micro_outside_MTOC_2 / ( double ) number_of_points_outside_MTOC_2 - sim_of_Cell::resting_distance  );
	    double b_b = abs( lenght_micro_outside_MTOC_2 / ( double ) ( number_of_points_outside_MTOC_2 - 1.0 ) - sim_of_Cell::resting_distance );

	   if( b_b <= a_a )
	   {
		this->resizeMicrotubule_minus_one();
	   }
           else
	   {
		this->resizeMicrotubule_same_number_of_beads();
	   }


       }
	//One bead can be added       
       else if( this->getTangent2( this->numberOfPoints - 2 ).norm() > 1.5 * sim_of_Cell::resting_distance  ) //1.2
       {
	    counter = counter + 1;
	    double lenght_micro_outside_MTOC_2 = this->get_lenght_of_microtubule_outside_MTOC();
	    unsigned int number_of_points_outside_MTOC_2 = this->getNumberOfPoints() - 2;
	    double a_a = abs( lenght_micro_outside_MTOC_2 / ( double ) number_of_points_outside_MTOC_2 - sim_of_Cell::resting_distance  );
	    double b_b = abs( lenght_micro_outside_MTOC_2 / ( double ) ( number_of_points_outside_MTOC_2 + 1.0 ) - sim_of_Cell::resting_distance );

	   if( a_a <= b_b )
	   {
		this->resizeMicrotubule_same_number_of_beads();	

	   }
           else
	   {
		this->resizeMicrotubule_plus_one();

	   }

       }
   }
   else
   {

   }

    this->set_lenght_of_tangents("void  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0()");
    double lenght_2 = this->get_lenght_of_microtubule_outside_MTOC();

    if( abs( lenght_2 - lenght_1 ) > 1e-9 )
    {	    
	cout<<"	 abs( lenght_2 - lenght_1 ) > 1e-9 "<<endl;	
	cout<<"lenght_2 = "<<lenght_2<<endl;
	cout<<"lenght_1 = "<<lenght_1<<endl;
	cout<<"void  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2_dynein_index_0()"<<endl;
	cout<<"   "<<endl;
	cout<<"ERROR_ID = 797979898"<<endl;
	throw("");
    }	




}









//Controls if the dyneins detach due to depolymerization
std::vector<Vector3d> Microtubule::control_motor_detachment_with_depolimerization()
{

    //If the dynein abscissa is longer than the microtubule - detachment
    std::vector< Vector3d > kicked_vector;
    std::vector< std::pair < Vector3d , double > > remaining_vector;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d ,double > pair_tmp = this->Dynein_motors_2[ pair_index ];
        double abscissa = std::get< 1 >( pair_tmp );
        Vector3d point = std::get< 0 >( pair_tmp );
	//This deals with imprecision emerging due to the resizing
        if( abs( abscissa - this->get_lenght_of_microtubule() ) < Dynein_real::step / 10.0 )
        {
            abscissa = abscissa - Dynein_real::step / 10.0;
        }

        std::pair < Vector3d ,double > pair_tmp_2( point , abscissa );


        if( abscissa > this->get_lenght_of_microtubule() )
        {
            kicked_vector.push_back( point );
        }
        else
        {
            remaining_vector.push_back( pair_tmp_2 );
        }
    }
    this->Dynein_motors_2 = remaining_vector;
    return kicked_vector;
}









//Returns the number of points attached to the microbule
unsigned int Microtubule::get_number_of_dynein_points_IS()
{
    return this->Dynein_motors_2.size();
}




//Controls the length of the microtubule and compares it with the distance between the MTOC and the depolymerization point 
std::vector<Vector3d> Microtubule::control_length_of_micro_IS_2()
{

    if( this->dydein_index != 20 )
    {
	cout<<"this->dydein_index != 20"<<endl;
	cout<<"ERROR_ID = 5468435615646"<<endl;
	throw("");
    }


    //Vector3d MTOC_center = ( this->getPoint( 0 ) + this->getPoint( 1 ) ) / 2.0;
    double distance = ( this->getPoint( 0 ) - this->IS_position_catching ).norm();

    for(  unsigned int bead_id = 2 ; bead_id < this->getNumberOfPoints() ;  bead_id ++ )
    {

	Vector3d bead_position = this->getPoint( bead_id );
	double distance_2 = ( bead_position - this->IS_position_catching ).norm();
	if( distance_2 > distance * 1.1 )
	{
	     std::vector<Vector3d> erased_vectors = this->get_dynein_points_in_IS_and_erase();

	     Vector3d posledni_tangent = this->get_last_Tangent();
	     Vector3d position = this->getPoint( this->getNumberOfPoints() - 2 );
	     position = position + posledni_tangent / posledni_tangent.norm() * sim_of_Cell::resting_distance;
	     this->setPoint( this->getNumberOfPoints() - 1 , position );
	     this->set_lenght_of_tangents("unsigned int Microtubule::get_number_of_dynein_points_IS()");

	     this->dydein_index	= 0;
	     return erased_vectors;
	}
    }

    std::vector<Vector3d> empty_vector;
    return empty_vector;
}



//Controls whether NaN exists in the variables of the microtubule 
void Microtubule::control_NAN(  string arg )
{

    for(  unsigned int bead_id = 0 ; bead_id < this->getNumberOfPoints() ;  bead_id ++ )
    {
    	//If the NaN is presented, the simulation is stopped
	Vector3d bead_position = this->getPoint( bead_id );
	if( std::isnan(  bead_position.norm()  )  )
	{
		cout<<arg<<endl;
		cout<<"bead_id = "<<bead_id<<endl;
		cout<<"bead_position.norm() = "<<bead_position.norm()<<endl;
		cout<<"Microtubule::control_NAN()"<<endl;
		cout<<" ERROR_ID = 3131331313 "<<endl;
		cout<<"cell_id = "<<this->cell_id<<endl;
		cout<<this->getCoordinates()<<endl;
		throw("");
	}
	bool answer = this->control_nan_Utilities( this->lenght_of_tangents_micro_0 );
	if( answer == true )
	{
		cout<<"bool answer = this->control_nan_Utilities( this->lenght_of_tangents_micro_0 )"<<endl;
		cout<<arg<<endl;
		cout<<this->lenght_of_tangents_micro_0<<endl;
		cout<<"cell_id = "<<this->cell_id<<endl;
		throw("");
	}
	
    }


}


//Returns all dyneins acting on a microtubule attached to cortical sliding
std::vector<Vector3d> Microtubule::get_catching_points_IS_sliding( )
{
    if( this->get_dynein_index() != 9 )
    {
        cout<<"this->get_dynein_index() != 9"<<endl;
        cout<<"std::vector< Vector3d  > Microtubule::get_Dynein_points_and_erase()"<<endl;
        cout<<"ERROR_ID = 6835461681364813"<<endl;
        throw("");
    }

    std::vector< Vector3d  > point_coordinates;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d anchor_position = std::get<0>( pair_tmp );
        point_coordinates.push_back( anchor_position );
    }
    return point_coordinates;
}












//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Detaches from dyneins and the index is set to zero
std::vector<Vector3d> Microtubule::detach_from_IS_2()
{
    std::vector<Vector3d> erased_vectors = this->get_dynein_points_in_IS_and_erase();
    this->dydein_index	= 0;
    if( this->get_last_Tangent().norm() < this->getRestDist() / 3 )	
    {

    }	
    return erased_vectors;

}


//Basic polymerization of the microtubule
void Microtubule::microtubule_growth_2()
{
	//tohle jen posoupne posledni bead mikrotubuly
	if( this->growth_index != 1 )
	{
		cout<<"this->growth_index != 1"<<endl;
		cout<<"void Microtubule::microtubule_growth_2()"<<endl;
		cout<<"ERROR_ID = 64568461"<<endl;
		throw("");
	}

	if( this->dydein_index == 0 )
	{
		this->resizeMicrotubule_with_different_tangent_lenghts();
		//Micro Growth
		if( this->numberOfPoints < 20 )
		{
			Vector3d last_tangent = this->get_last_Tangent();
			Vector3d last_point = this->getPoint( this->numberOfPoints - 2 );
			Vector3d new_last_bead = last_point + last_tangent / last_tangent.norm() *  ( last_tangent.norm() + dynamic_instability::polymerization_constant * sim_of_Cell::time_Step );
			for( unsigned int index = 0 ; index < 3 ; index ++ )
			{
				this->coordinates( 3 * ( this->numberOfPoints - 1 ) + index , 0 ) = new_last_bead( index );
			}
		}

	}
}


//This is the growth of the microtubule due to the dynamic instability
void Microtubule::microtubule_growth_MT_0_9()
{
	//tohle jen posoupne posledni bead mikrotubuly
	if( this->growth_index != 1 )
	{
		cout<<"this->growth_index != 1"<<endl;
		cout<<"void Microtubule::microtubule_growth_MT_0_9()"<<endl;
		cout<<"ERROR_ID = 64568461"<<endl;
		throw("");
	}
	//Only the microtubule with a free end can grow
	if( ( this->dydein_index != 0 ) && ( this->dydein_index != 9 ) )
	{
		cout<<" ( this->dydein_index != 0 ) && ( this->dydein_index != 9 ) "<<endl;
		cout<<"void Microtubule::microtubule_growth_MT_0_9()"<<endl;
		cout<<"ERROR_ID = 4646466446446"<<endl;
		throw("");
	}


	double last_tangent_0 = this->get_last_Tangent().norm();
	//grow of the microtubule
	if( this->dydein_index == 0 )
	{
		this->resizeMicrotubule_with_different_tangent_lenghts();
		Vector3d last_tangent = this->get_last_Tangent();
		Vector3d last_point = this->getPoint( this->numberOfPoints - 2 );
		Vector3d new_last_bead = last_point + last_tangent / last_tangent.norm() *  ( last_tangent.norm() + dynamic_instability::polymerization_constant * sim_of_Cell::time_Step );
		for( unsigned int index = 0 ; index < 3 ; index ++ )
		{
			this->coordinates( 3 * ( this->numberOfPoints - 1 ) + index , 0 ) = new_last_bead( index );
		}


	}
	else if( this->dydein_index == 9 )
	{
		this->resizeMicrotubule_with_different_tangent_lenghts();
		Vector3d last_tangent = this->get_last_Tangent();
		Vector3d last_point = this->getPoint( this->numberOfPoints - 2 );
		Vector3d new_last_bead = last_point + last_tangent / last_tangent.norm() *  ( last_tangent.norm() + dynamic_instability::polymerization_constant * sim_of_Cell::time_Step );
		for( unsigned int index = 0 ; index < 3 ; index ++ )
		{
			this->coordinates( 3 * ( this->numberOfPoints - 1 ) + index , 0 ) = new_last_bead( index );
		}


	}
	this->set_lenght_of_tangents("void Microtubule::microtubule_growth_MT_0_9()");
	double last_tangent_1 = this->get_last_Tangent().norm();


}




//Shrinking of the microtubule with a free end
void Microtubule::microtubule_shrinkage_MT_0_9()
{
	//tohle jen posoupne posledni bead mikrotubuly
	if( this->growth_index != 2 )
	{
		cout<<"this->growth_index != 2"<<endl;
		cout<<"void Microtubule::microtubule_shrinkage_MT_0_9()"<<endl;
		cout<<"ERROR_ID = 98796844964"<<endl;
		throw("");
	}
	if( ( this->dydein_index != 0 ) && ( this->dydein_index != 9 ) )
	{
		cout<<" ( this->dydein_index != 0 ) && ( this->dydein_index != 9 ) "<<endl;
		cout<<"void Microtubule::microtubule_shrinkage_MT_0_9()"<<endl;
		cout<<"ERROR_ID = 161565615661"<<endl;
		throw("");
	}

	double last_tangent_0 = this->get_last_Tangent().norm();
	//Shrinkage of the microtubule
	if( this->dydein_index == 0 )
	{
		double aaa = this->get_last_Tangent().norm();
		this->resizeMicrotubule_with_different_tangent_lenghts();
		Vector3d last_tangent = this->get_last_Tangent();
		Vector3d last_point = this->getPoint( this->numberOfPoints - 2 );
		Vector3d new_last_bead = last_point + last_tangent / last_tangent.norm() *  ( last_tangent.norm() - dynamic_instability::shrinking_constant * sim_of_Cell::time_Step );
		for( unsigned int index = 0 ; index < 3 ; index ++ )
		{
			this->coordinates( 3 * ( this->numberOfPoints - 1 ) + index , 0 ) = new_last_bead( index );
		}
		double bbb = this->get_last_Tangent().norm();

	}
	if( this->dydein_index == 9 )
	{
		this->resizeMicrotubule_with_different_tangent_lenghts();
		Vector3d last_tangent = this->get_last_Tangent();
		Vector3d last_point = this->getPoint( this->numberOfPoints - 2 );
		Vector3d new_last_bead = last_point + last_tangent / last_tangent.norm() *  ( last_tangent.norm() - dynamic_instability::shrinking_constant * sim_of_Cell::time_Step );
		for( unsigned int index = 0 ; index < 3 ; index ++ )
		{
			this->coordinates( 3 * ( this->numberOfPoints - 1 ) + index , 0 ) = new_last_bead( index );
		}

	}
	this->set_lenght_of_tangents("void Microtubule::microtubule_shrinkage_MT_0_9()");

}



















//Controls whether the microtubule was overextended due to the repeated attachment to capture-shrinkage dynein where the tip is placed on the membrane
bool Microtubule::control_of_captured_MT_lenght_against_prolongation()
{
	//return false;
	double chosen_prolongation_treshold = 1.0 * sim_of_Cell::chosen_prolongation_tr;
	if( this->get_lenght_of_microtubule_outside_MTOC(  ) > this->get_lenght_of_micro_0_outside_MTOC(  ) + chosen_prolongation_treshold ) 
	{
		return true;
	}
	else
	{
		return false;
	}
}













//Shortens the microtubule if it became overextended due to repeated attachment
void Microtubule::shorten_MT_due_to_numerical_prolongation_during_attachment()
{

	double lenght_1 = this->get_lenght_of_microtubule_outside_MTOC(); 	
	double shortening = this->get_lenght_of_microtubule_outside_MTOC(  ) - get_lenght_of_micro_0_outside_MTOC(  );

	if( shortening < 0 )
	{
		cout<<"  shortening < 0  "<<endl;
		cout<<"void Microtubule::shorten_MT_due_to_numerical_prolongation_during_attachment()"<<endl;
		cout<<"ERROR_ID = 111611"<<endl;
		throw("");

	}


	if( this->get_last_Tangent( ).norm() <  shortening )
	{	
		double shortening_tmp = shortening - this->get_last_Tangent( ).norm();
		this->Substract_one_Bead_simple(  );


		Vector3d point_before_last = this->getPoint( this->getNumberOfPoints() - 2 );
		Vector3d last_Tangent = this->get_last_Tangent( );
		Vector3d last_Point = point_before_last + last_Tangent / last_Tangent.norm() * ( last_Tangent.norm() -  shortening_tmp );

		for( unsigned int index_tmp = 0 ; index_tmp < 3 ; index_tmp ++ )
		{
			this->coordinates( 3 * ( this->getNumberOfPoints() - 1 ) + index_tmp ,  0 ) = last_Point( index_tmp );
		}
                //This justmakes sure that the last segment is not too short
                //Controls will adjust the segments anyway
		if( this->get_last_Tangent( ).norm() < 0.33 * this->getRestDist() ) 
		{
			if( this->getNumberOfPoints() > 3 )
			{
				unsigned int number_of_points = this->getNumberOfPoints();
				Vector3d last_tangent = this->get_last_Tangent( );
				Vector3d previous_to_last_Tangent = this->getTangent2( number_of_points - 3 );
				double lenght_of_new_last_tangent = last_tangent.norm() + previous_to_last_Tangent.norm();
				Vector3d point_two_before_last = this->getPoint( this->getNumberOfPoints() - 3 );
				Vector3d new_last_point = point_two_before_last + previous_to_last_Tangent / previous_to_last_Tangent.norm() * ( lenght_of_new_last_tangent ); 


				this->Substract_one_Bead_simple(  );
				for( unsigned int index_tmp = 0 ; index_tmp < 3 ; index_tmp ++ )
				{
					this->coordinates( 3 * ( this->getNumberOfPoints() - 1 ) + index_tmp ,  0 ) = new_last_point( index_tmp );
				}
			 }
		}
	}
	else
	{
		Vector3d point_before_last = this->getPoint( this->getNumberOfPoints() - 2 );
		Vector3d last_Tangent = this->get_last_Tangent( );
		Vector3d last_Point = point_before_last + last_Tangent / last_Tangent.norm() * ( last_Tangent.norm() -  shortening );
		for( unsigned int index_tmp = 0 ; index_tmp < 3 ; index_tmp ++ )
		{
			this->coordinates( 3 * ( this->getNumberOfPoints() - 1 ) + index_tmp ,  0 ) = last_Point( index_tmp );
		}

		if( this->get_last_Tangent( ).norm() < 0.33 * this->getRestDist() ) 
		{
			if( this->getNumberOfPoints() > 3 )
			{
				unsigned int number_of_points = this->getNumberOfPoints();
				Vector3d last_tangent = this->get_last_Tangent( );
				Vector3d previous_to_last_Tangent = this->getTangent2( number_of_points - 3 );
				double lenght_of_new_last_tangent = last_tangent.norm() + previous_to_last_Tangent.norm();
				Vector3d point_two_before_last = this->getPoint( this->getNumberOfPoints() - 3 );
				Vector3d new_last_point = point_two_before_last + previous_to_last_Tangent / previous_to_last_Tangent.norm() * ( lenght_of_new_last_tangent ); 


				this->Substract_one_Bead_simple(  );
				for( unsigned int index_tmp = 0 ; index_tmp < 3 ; index_tmp ++ )
				{
					this->coordinates( 3 * ( this->getNumberOfPoints() - 1 ) + index_tmp ,  0 ) = new_last_point( index_tmp );
				}
			 }
		}

	}

	this->set_lenght_of_tangents( "void Microtubule::shorten_MT_due_to_numerical_prolongation_during_attachment()" );	

        //This is just a control of the function
	double lenght_2 = this->get_lenght_of_microtubule_outside_MTOC(); 
	double difference = lenght_2 - lenght_1;
	if( abs( abs( difference ) - shortening  ) >  1e-9  )
	{
		double difference = lenght_2 - lenght_1;
		cout<<"difference = "<<difference<<endl;
		cout<<"shortening = "<<shortening<<endl;
		cout<<"abs( difference - shortening ) = "<<abs( difference - shortening )<<endl;
		cout<<"void Microtubule::shorten_MT_due_to_numerical_prolongation_during_attachment()"<<endl;
		cout<<"ERROR_ID = 74684464646"<<endl;
		throw("");
	} 
	






}

































//Resizes growing, unattached microtubule
void  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0()
{
    //Carefull, do not insert grow!
    //It is used only when the micro is growing
    //Therefore, it has to be resized before growing

    if( ( this->get_dynein_index() != 0 ) && ( this->get_dynein_index() != 9 ) )
    {
        cout<<"Microtubule::control_and_resizing_of_depolimerizing_microtubule()"<<endl;
        cout<<" ( this->get_dynein_index() != 0 ) "<<endl;
        unsigned int ERROR_ID = 614611;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
	cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        throw ERROR_ID;
    }
    double lenght_1 = this->get_lenght_of_microtubule_outside_MTOC(); 	

    if( this->numberOfPoints <= 2  )
    {
        cout<<"this->numberOfPoints <= 2"<<endl;
        cout<<"void  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0()"<<endl;
	cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
	cout<<"ERROR_ID = 9798798779797979779"<<endl;
	throw("");

    }


    if( this->get_last_Tangent().norm()  > 1.2 * sim_of_Cell::resting_distance )
    {
	    double lenght_micro_outside_MTOC_2 = this->get_lenght_of_microtubule_outside_MTOC();
	    unsigned int number_of_points_outside_MTOC_2 = this->getNumberOfPoints() - 2;
	    double a_a = abs( lenght_micro_outside_MTOC_2 / ( double ) number_of_points_outside_MTOC_2 - sim_of_Cell::resting_distance  );
	    double b_b = abs( lenght_micro_outside_MTOC_2 / ( double ) ( number_of_points_outside_MTOC_2 + 1.0 ) - sim_of_Cell::resting_distance );

	   //This decides whether an additional bead is created
	   if( a_a <= b_b )
	   {
		this->resizeMicrotubule_same_number_of_beads();
	   }
           else
	   {
		this->resizeMicrotubule_plus_one();
	   }
	   double friction = this->calculateEffectiveFriction_Howard( );
	   this->set_effective_friction( friction );		
   }
	
    double lenght_2 = this->get_lenght_of_microtubule_outside_MTOC();
    if( abs( lenght_2 - lenght_1 ) > 1e-9 )
    {	    
	cout<<"	 abs( lenght_2 - lenght_1 ) > 1e-9 "<<endl;	
	cout<<"lenght_2 = "<<lenght_2<<endl;
	cout<<"lenght_1 = "<<lenght_1<<endl;
	cout<<"void  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0()"<<endl;
	cout<<"   "<<endl;
	cout<<"ERROR_ID = 313513351313131313"<<endl;
	throw("");
    }	

    this->set_lenght_of_tangents("void  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_dynein_index_0()");
}


//Controls a shrinking microtubule and resizes it
void  Microtubule::control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9()
{
    //This function cannot contain resizing
    //Only use when the microtubule is depolymerized

    if( ( this->get_dynein_index() != 0 ) && ( this->get_dynein_index() != 9 )  )
    {

        cout<<"control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9"<<endl;
        cout<<" ( this->get_dynein_index() != 0 ) && ( this->get_dynein_index() != 9 )  "<<endl;
        unsigned int ERROR_ID = 351561456;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }

    if( this->numberOfPoints <= 2  )
    {
        cout<<"this->numberOfPoints <= 2"<<endl;
        cout<<"void  Microtubule::control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9()"<<endl;
	cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
	cout<<"ERROR_ID = 22113131313131313131"<<endl;
	throw("");
    }
    double lenght_1 = this->get_lenght_of_microtubule_outside_MTOC();


    if( this->numberOfPoints > 3 )
    {

       if( this->getTangent2( this->numberOfPoints - 2 ).norm() < 0.8 * sim_of_Cell::resting_distance  )
       {
	    double lenght_micro_outside_MTOC_2 = this->get_lenght_of_microtubule_outside_MTOC();
	    unsigned int number_of_points_outside_MTOC_2 = this->getNumberOfPoints() - 2;
	    double a_a = abs( lenght_micro_outside_MTOC_2 / ( double ) number_of_points_outside_MTOC_2 - sim_of_Cell::resting_distance  );
	    double b_b = abs( lenght_micro_outside_MTOC_2 / ( double ) ( number_of_points_outside_MTOC_2 - 1.0 ) - sim_of_Cell::resting_distance );
	   //One bead is removed
	   if( b_b <= a_a )
	   {
		this->resizeMicrotubule_minus_one();
	   }
           else
	   {
		this->resizeMicrotubule_same_number_of_beads();
	   }
	   double friction = this->calculateEffectiveFriction_Howard( );
	   this->set_effective_friction( friction );

       }

    }
    this->set_lenght_of_tangents("void  Microtubule::control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9()");

    //The control of the function
    double lenght_2 = this->get_lenght_of_microtubule_outside_MTOC();
    if( abs( lenght_2 - lenght_1 ) > 1e-9 )
    {	    
	cout<<"	 abs( lenght_2 - lenght_1 ) > 1e-9 "<<endl;	
	cout<<"lenght_2 = "<<lenght_2<<endl;
	cout<<"lenght_1 = "<<lenght_1<<endl;
	cout<<"void  Microtubule::control_and_resizing_and_motor_adjustment_SHRINKING_microtubule_dynein_index_0_9()"<<endl;
	cout<<"   "<<endl;
	cout<<"ERROR_ID = 616461561651651355661566616"<<endl;
	throw("");
    }	

}






//It resizes the microtubule so that all segments have the same size and the number of beads stays the same
void Microtubule::resizeMicrotubule_same_number_of_beads()
{

    if( ( this->get_dynein_index() != 0 ) && ( this->get_dynein_index() != 9 ) )
    {
        cout<<"Microtubule::resizeMicrotubule_same_number_of_beads()"<<endl;
        cout<<" ( this->get_dynein_index() != 0 ) "<<endl;
        unsigned int ERROR_ID = 322131313;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
	cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        throw ERROR_ID;
    }
    double lenght_1 = this->get_lenght_of_microtubule_outside_MTOC(); 	

    if( this->numberOfPoints <= 2  )
    {
        cout<<"this->numberOfPoints <= 2"<<endl;
        cout<<"void  Microtubule::resizeMicrotubule_same_number_of_beads()"<<endl;
	cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
	cout<<"ERROR_ID = 456464564564465"<<endl;
	throw("");

    }

	//Saves the orietation of the segments
	std::vector< Vector3d > orientations;
    	for( unsigned int point_id = 1 ; point_id < this->getNumberOfPoints() - 1 ; point_id ++ )
    	{
		orientations.push_back( this->getTangent2( point_id ) / this->getTangent2( point_id ).norm() );
    	}
    	double lenght_micro_outside_MTOC = this->get_lenght_of_microtubule_outside_MTOC();

    	unsigned int number_of_points_outside_MTOC = this->getNumberOfPoints() - 2;
    	//New length of the segment
	if( this->getNumberOfPoints() > 3 )
	{
    		this->restDistancePoints = lenght_micro_outside_MTOC / ( double ) ( number_of_points_outside_MTOC );
	}
	else
	{
		this->restDistancePoints = lenght_micro_outside_MTOC;
    	}

       	MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    	for( unsigned int counter = 0 ; counter < 6 ; counter ++ )
    	{
		tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
    	}

   	for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
   	{
		Vector3d orient = this->restDistancePoints * orientations[ i - 2 ];
		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
			tmp_new( 3 * i + j , 0 ) = tmp_new( 3 * ( i - 1 ) + j , 0 ) + orient( j );
		}
   	}
    	this->coordinates = tmp_new;
}



//It resizes the microtubule so that all segments have the same size, number of beads is increased
void Microtubule::resizeMicrotubule_plus_one()
{

    if( ( this->get_dynein_index() != 0 ) && ( this->get_dynein_index() != 9 ) )
    {
        cout<<"Microtubule::resizeMicrotubule_plus_one()"<<endl;
        cout<<" ( this->get_dynein_index() != 0 ) "<<endl;
        unsigned int ERROR_ID = 21144414;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
	cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        throw ERROR_ID;
    }
    double lenght_1 = this->get_lenght_of_microtubule_outside_MTOC(); 	

    if( this->numberOfPoints <= 2  )
    {
        cout<<"this->numberOfPoints <= 2"<<endl;
        cout<<"void  Microtubule::resizeMicrotubule_plus_one()"<<endl;
	cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
	cout<<"ERROR_ID = 17178178871781"<<endl;
	throw("");

    }


	//Length of the microtubule outside the MTOC
	double lenght_micro_outside_MTOC = this->get_lenght_of_microtubule_outside_MTOC();
 	//Orientations are saved
	std::vector< Vector3d > orientations;
    	for( unsigned int point_id = 1 ; point_id < this->getNumberOfPoints() - 1 ; point_id ++ )
    	{
		orientations.push_back( this->getTangent2( point_id ) / this->getTangent2( point_id ).norm() );
    	}
	//Create one additional orientation
	Vector3d last_orienation = this->get_last_Tangent( ) / this->get_last_Tangent( ).norm();
	orientations.push_back( last_orienation );


    	//Change the number of beads
   	this->numberOfPoints = this->numberOfPoints + 1;
    	unsigned int number_of_points_outside_MTOC = this->getNumberOfPoints() - 2;
    	//It is already changed for shorter microtubule
    	this->restDistancePoints = lenght_micro_outside_MTOC / ( double ) ( number_of_points_outside_MTOC );

	//First two beads are unchanged
       	MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    	for( unsigned int counter = 0 ; counter < 6 ; counter ++ )
    	{
		tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
    	}

   	for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
   	{
		Vector3d orient = this->restDistancePoints * orientations[ i - 2 ];
		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
			tmp_new( 3 * i + j , 0 ) = tmp_new( 3 * ( i - 1 ) + j , 0 ) + orient( j );
		}
   	}

   	this->coordinates = tmp_new;
}
//Resizes the microtubule so that all segments have the same size, number of beads is reduced
void Microtubule::resizeMicrotubule_minus_one()
{
    if( ( this->get_dynein_index() != 0 ) && ( this->get_dynein_index() != 9 ) )
    {
        cout<<"Microtubule::resizeMicrotubule_minus_one()"<<endl;
        cout<<" ( this->get_dynein_index() != 0 ) "<<endl;
        unsigned int ERROR_ID = 363636336;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
	cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        throw ERROR_ID;
    }
    double lenght_1 = this->get_lenght_of_microtubule_outside_MTOC(); 	

    if( this->numberOfPoints <= 2  )
    {
        cout<<"this->numberOfPoints <= 2"<<endl;
        cout<<"void  Microtubule::resizeMicrotubule_minus_one()"<<endl;
	cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
	cout<<"ERROR_ID = 36363626262"<<endl;
	throw("");
    }
		std::vector< Vector3d > orientations;
	    	for( unsigned int point_id = 1 ; point_id < this->getNumberOfPoints() - 1 ; point_id ++ )
	    	{
			orientations.push_back( this->getTangent2( point_id ) / this->getTangent2( point_id ).norm() );
	    	}
	    	double lenght_micro_outside_MTOC = this->get_lenght_of_microtubule_outside_MTOC();

	    	//Change of the number of beads
	   	this->numberOfPoints = this->numberOfPoints - 1;
	    	unsigned int number_of_points_outside_MTOC = this->getNumberOfPoints() - 2;
	    	//New length of the segment
    		this->restDistancePoints = lenght_micro_outside_MTOC / ( double ) ( number_of_points_outside_MTOC );

            	MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    	    	for( unsigned int counter = 0 ; counter < 6 ; counter ++ )
    	    	{
			tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
    	    	}

    	   	for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
    	   	{
			Vector3d orient = this->restDistancePoints * orientations[ i - 2 ];
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				tmp_new( 3 * i + j , 0 ) = tmp_new( 3 * ( i - 1 ) + j , 0 ) + orient( j );
			}
    	   	}

  	  	this->coordinates = tmp_new;
}





//Controls whether NaN exist in the matrix
bool Microtubule::control_nan_Utilities( MatrixXd tmp  )
{
	if( tmp.cols() != 1 )
	{
		cout<<" tmp.cols() != 1 "<<endl;
		cout<<"void control_nan_Utilities( MatrixXd tmp )"<<endl;
		throw("");
	}	

	for( unsigned int point_tmp = 0 ; point_tmp < tmp.rows() ; point_tmp ++ )
	{
		if(  std::isnan( tmp(  point_tmp ) ) == true  )
		{
			cout<<tmp<<endl;
			cout<<"cell_id = "<<this->cell_id<<endl;
			return true;
		} 
	}
	return false;

}









Microtubule::~Microtubule()
{
	// TODO Auto-generated destructor stub
}



//////////////////////////////////////////////////////////////////////////////Unused constructors

Microtubule::Microtubule( Vector3d MTOC , unsigned int ID )
{
	this->cell_id = 0;
	//simple default constructor - direct microtubule just with id
	//on x axis - just starts in MTOC point
	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
    	this->MTOC_point = 1;
	this->polygon = 10000;
	//radius has no meaning - the cell do not restrict Microtubule
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
	this->effective_friction = this->calculateEffectiveFriction_Howard();
	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 0.0 , 0.0 , 0.0 );

	if( this->numberOfPoints < 1  )
	{
		cout<<" this->numberOfPoints < 1 in Microtubule::Microtubule( Vector3d MTOC , unsigned int ID )"<<endl;
		throw("");
	}

	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );

	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		this->coordinates( 3 * i , 0 ) = MTOC( 0 ) +  i * this->restDistancePoints;
		this->coordinates( 3 * i + 1 , 0 ) = MTOC( 1 ) + 0.0;
		this->coordinates( 3 * i + 2 , 0 ) = MTOC( 2 ) + 0.0;
	}
	//this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
    this->set_growth_index( 0 );
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	
}



Microtubule::Microtubule( Vector3d MTOC , Vector3d first_Point , unsigned int ID )
{
//	cout<<"ID = "<<ID<<endl;
//	cout<<first_Point<<endl;
	this->cell_id = 0;
	if( ( MTOC - first_Point ).norm() == 0 )
	{
		cout<<"( MTOC - first_Point ).norm() == 0 in Microtubule::Microtubule( Vector3d MTOC , Vector3d first_Point , unsigned int ID )"<<endl;
		throw("");
	}


	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
	this->polygon = 10000;
    this->MTOC_point = 1;
	//radius has no meaning - the cell do not restrict Microtubule
	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 0.0 , 0.0 , 0.0 );
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
	this->effective_friction = this->calculateEffectiveFriction_Howard();
        //this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );

	if( this->numberOfPoints < 1  )
	{
		cout<<" this->numberOfPoints < 1 in Microtubule::Microtubule( Vector3d MTOC , Vector3d first_Point , unsigned int ID )"<<endl;
		throw("");
	}



	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );

	//direction from MTOC to first_point
	//this is direction of microtubule in x - y plane - does not change for entire microtubule
	Vector3d direction( ( first_Point[ 0 ] - MTOC[ 0 ]  ) , ( first_Point[ 1 ]  - MTOC[ 1 ] )  , 0 );
	direction = direction / direction.norm();

	//axis that I will use to rotate the vectors
	Vector3d axis( direction( 1 ) , - direction( 0 )  , 0 );

	//constructing Quaterion - object to rotate vectors
	//the same Quaterion will be used for all rotation
	Quaternion<double> q;
	q = AngleAxis<double>( sim_of_Cell::angle , axis );



	//coordinates of first point:: equals first_Point
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );


	if( this->numberOfPoints > 1 )
	{
		//first tangent - between first_Point and second bead of microtubule
		//REMARK!!! microtubule is not influenced by distance of MTOC to first point
		Vector3d tangent( this->restDistancePoints * direction( 0 ) , this->restDistancePoints * direction( 1 ) , 0.0 );

		this->coordinates( 3 , 0 ) = first_Point( 0 ) + tangent( 0 );
		this->coordinates( 4 , 0 ) = first_Point( 1 ) + tangent( 1 );
		this->coordinates( 5 , 0 ) = first_Point( 2 ) + tangent( 2 );


		//Others tangents are created here- there magnitude is still the same
		//Only difference is orientation
		//REMEMBER - orientation in xy plane remains the same: x / y = constant -
		//that is because "spider-like" structure of microtubule cytoskeleton - see movies
		for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
		{
			tangent = q * tangent;
			this->coordinates( 3 * i , 0 ) = this->coordinates( 3 * ( i - 1 ) , 0 ) + tangent( 0 );
			this->coordinates( 3 * i + 1 , 0 ) = this->coordinates( 3 * ( i - 1 ) + 1 , 0 ) + tangent( 1 );
			this->coordinates( 3 * i + 2 , 0 ) = this->coordinates( 3 * ( i - 1 ) + 2 , 0 ) + tangent( 2 );
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
    this->set_growth_index( 0 );
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	
}



Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID , unsigned int poly )
{
	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points
	//TMP has no meaning and is used only to overload function
	this->cell_id = 0;
	if( ( second_Point - first_Point ).norm() == 0 )
	{
		cout<<"( second_Point - first_Point ).norm() == 0 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID , 666 )"<<endl;
		throw("");
	}

	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->MTOC_point = 1;
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
	this->effective_friction = this->calculateEffectiveFriction_Howard();
	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );
    //this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );

	if( abs( ( second_Point - first_Point ).norm() / this->restDistancePoints - 1.0 ) > 1e-8)
	{
		cout<<"abs( ( second_Point - first_Point ).norm() / this->restDistancePoints ) > 1e-8 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID , 666 )"<<endl;
		throw("");
	}


	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );

	//direction from MTCO to first_point
	//this is direction of microtubule in x - y plane - does not change for entire microtubule
	Vector3d direction( second_Point[ 0 ] - first_Point[ 0 ] , second_Point[ 1 ] - first_Point[ 1 ]  , second_Point[ 2 ] - first_Point[ 2 ] );
	direction = direction / direction.norm();


	//coordinates of first point:: equals first_Point
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

	this->coordinates( 3 , 0 ) = second_Point( 0 );
	this->coordinates( 4 , 0 ) = second_Point( 1 );
	this->coordinates( 5 , 0 ) = second_Point( 2 );


	if( this->numberOfPoints > 2 )
	{
		//first tangent - between first_Point and second bead of microtubule
		//REMARK!!! microtubule is not influenced by distance of MTOC to first point
		Vector3d tangent( this->restDistancePoints * direction( 0 ) , this->restDistancePoints * direction( 1 ) , 0.0 );


		//Others tangents are created here- there magnitude is still the same
		//Only difference is orientation
		//REMEMBER - orientation in xy plane remains the same: x / y = constant -
		//that is because "spider-like" structure of microtubule cytoskeleton - see movies
		for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
		{
			this->coordinates( 3 * i , 0 ) = this->coordinates( 3 * ( i - 1 ) , 0 ) + tangent( 0 );
			this->coordinates( 3 * i + 1 , 0 ) = this->coordinates( 3 * ( i - 1 ) + 1 , 0 ) + tangent( 1 );
			this->coordinates( 3 * i + 2 , 0 ) = this->coordinates( 3 * ( i - 1 ) + 2 , 0 ) + tangent( 2 );
		}
	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
    this->set_growth_index( 0 );
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	
}


Microtubule::Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , double a_axis , double b_axis )
{
	this->cell_id = 0;
	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points
	//TMP has no meaning and is used only to overload function

	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->MTOC_point = 1;
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
	this->effective_friction = this->calculateEffectiveFriction_Howard();
    //this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );

	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );

	//confirm that first bead lies in the ellipsoid - if not throw error
	bool index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );




	if( index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"index = "<<index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=6789468416351684"<<endl;
		throw("");
	}


	//direction from MTOC to first_point
	//this is direction of microtubule in x - y plane - does not change for entire microtubule
	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	//axis that I will use to rotate the vectors
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	//coordinates of first point:: equals first_Point
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );


	//last point and orientation will be changed as we go from
	Vector3d last_point = first_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 1 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{
			//I aplly quaternion to rotate
			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;

			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
    this->set_growth_index( 0 );
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	

}






Microtubule::Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , double a_axis , double b_axis , unsigned int number_of_points  )
{
	this->cell_id = 0;
	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points
	//TMP has no meaning and is used only to overload function

	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = number_of_points;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->MTOC_point = 1;
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
	this->effective_friction = this->calculateEffectiveFriction_Howard();
    //this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );

	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}


	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );

	//confirm that first bead lies in the ellipsoid - if not throw error
	bool index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );




	if( index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"index = "<<index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=689451765"<<endl;
		throw("");
	}


	//direction from MTOC to first_point
	//this is direction of microtubule in x - y plane - does not change for entire microtubule
	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	//axis that I will use to rotate the vectors
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );

	//here I set the axis that will be used to get the axis for rotation
//	Quaternion<double> q_tmp;
//	q_tmp = AngleAxis<double>( rotation_angle , axis_rotation );

//	Vector3d axis_rotation = first_Point.cross( orientation );
//	axis_rotation = axis_rotation / axis_rotation.norm();


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	//coordinates of first point:: equals first_Point
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );


	//last point and orientation will be changed as we go from
	Vector3d last_point = first_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 1 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{
			//I aplly quaternion to rotate
			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;

			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
    this->set_growth_index( 0 );
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	
}


Microtubule::Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int number_of_points  )
{
	this->cell_id = 0;
    	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points
	//TMP has no meaning and is used only to overload function

	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = number_of_points;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->side = side_arg;
    this->MTOC_point = MTOC_point_arg;
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
this->effective_friction = this->calculateEffectiveFriction_Howard();
    //this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );
	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );
    double a_axis = Cell_parametres::A_AXIS;
    double b_axis = Cell_parametres::B_AXIS;
    //confirm that first bead lies in the ellipsoid - if not throw error
	bool bool_index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );

	if( bool_index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"bool_index = "<<bool_index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=46461343154127"<<endl;
		throw("");
	}


	//direction from MTOC to first_point
	//this is direction of microtubule in x - y plane - does not change for entire microtubule
	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	//axis that I will use to rotate the vectors
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );

	//here I set the axis that will be used to get the axis for rotation
//	Quaternion<double> q_tmp;
//	q_tmp = AngleAxis<double>( rotation_angle , axis_rotation );

//	Vector3d axis_rotation = first_Point.cross( orientation );
//	axis_rotation = axis_rotation / axis_rotation.norm();


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	//coordinates of first point:: equals first_Point
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );


	//last point and orientation will be changed as we go from
	Vector3d last_point = first_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 1 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{
			//I aplly quaternion to rotate
			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;

			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
    this->set_growth_index( 0 );
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	

}

Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points , bool rovna_mikrotubula )
{
	this->cell_id = 0;
    if( MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 689645684316464"<<endl;
        throw("");
    }
    if( second_MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 8987625148688"<<endl;
        throw("");
    }

    	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points
	//TMP has no meaning and is used only to overload function
    //this is the microtubule with first segment shorter
	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = number_of_points;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->side = side_arg;
    this->MTOC_point = MTOC_point_arg;
    this->second_MTOC_point = second_MTOC_point_arg;
    this->dydein_index = 0.0;
    this->effective_friction = this->calculateEffectiveFriction_Howard();
	//cout<<"this->calculateEffectiveFriction( ) = "<<this->calculateEffectiveFriction( )<<endl;
	//cout<<"this->calculateEffectiveFriction_Howard( ) = "<<this->calculateEffectiveFriction_Howard( )<<endl;
    //this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );
    this->restDistancePoints_first = 2.0 * MTOCparam::MTOC_radius;

    if( abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10  )
    {
        cout<<"( first_Point - second_Point ).norm() = "<<( first_Point - second_Point ).norm()<<endl;
        cout<<"this->restDistancePoints_first = "<<this->restDistancePoints_first<<endl;
        cout<<"abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10"<<endl;
        cout<<"ERROR_ID = 154355454435"<<endl;
        throw("");
    }



	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );
    double a_axis = Cell_parametres::A_AXIS;
    double b_axis = Cell_parametres::B_AXIS;
    //confirm that first bead lies in the ellipsoid - if not throw error
	bool bool_index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );

	if( bool_index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"bool_index = "<<bool_index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=46461343154127"<<endl;
		throw("");
	}


	//direction from MTOC to first_point
	//this is direction of microtubule in x - y plane - does not change for entire microtubule
	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	//axis that I will use to rotate the vectors
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        	this->coordinates( 3 + dimension , 0 ) = second_Point( dimension );
    }




	//last point and orientation will be changed as we go from
	Vector3d last_point = second_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 2 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{

			tmp_point = last_point + orientation * this->restDistancePoints;

			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }

    this->set_bending_matrix( );
    this->set_growth_index( 0 );
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	
}





Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points , double first_bead_distance )
{

	this->cell_id = 0;

    if( MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 689645684316464"<<endl;
        throw("");
    }
    if( second_MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 8987625148688"<<endl;
        throw("");
    }

    	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points
	//TMP has no meaning and is used only to overload function
    //this is the microtubule with first segment shorter
	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  sim_of_Cell::k_bending_analytical / restDistancePoints;
	this->numberOfPoints = number_of_points;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->side = side_arg;
    this->MTOC_point = MTOC_point_arg;
    this->second_MTOC_point = second_MTOC_point_arg;
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
	//cout<<"this->calculateEffectiveFriction( ) = "<<this->calculateEffectiveFriction( )<<endl;
	//cout<<"this->calculateEffectiveFriction_Howard( ) = "<<this->calculateEffectiveFriction_Howard( )<<endl;
	this->effective_friction = this->calculateEffectiveFriction_Howard();
    //this->Dynein_on_surface[ 0 ] = Vector3d( 0.0 , 0.0 , 0.0 );
    this->restDistancePoints_first = first_bead_distance;


    if( abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10  )
    {
        cout<<"( first_Point - second_Point ).norm() = "<<( first_Point - second_Point ).norm()<<endl;
        cout<<"this->restDistancePoints_first = "<<this->restDistancePoints_first<<endl;
        cout<<"abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10"<<endl;
        cout<<"ERROR_ID = 61543554545"<<endl;
        throw("");
    }



	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );
    double a_axis = Cell_parametres::A_AXIS;
    double b_axis = Cell_parametres::B_AXIS;
    //confirm that first bead lies in the ellipsoid - if not throw error
	bool bool_index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );

	if( bool_index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"bool_index = "<<bool_index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=46461343154127"<<endl;
		throw("");
	}


	//direction from MTOC to first_point
	//this is direction of microtubule in x - y plane - does not change for entire microtubule
	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	//axis that I will use to rotate the vectors
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        	this->coordinates( 3 + dimension , 0 ) = second_Point( dimension );
    }




	//last point and orientation will be changed as we go from
	Vector3d last_point = second_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 2 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}

			last_point = tmp_point;
		}
		else
		{
			//I aplly quaternion to rotate
			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;


			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}


    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    this->set_lenght_of_tangents("constructor");

    this->set_bending_matrix( );
    this->set_growth_index( 0 );
    this->set_lenght_of_micro_0_outside_MTOC( this->get_lenght_of_microtubule_outside_MTOC() );	
    this->set_lenght_of_tangents_micro_0();	
}





