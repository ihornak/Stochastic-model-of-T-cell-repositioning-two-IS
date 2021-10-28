/*
 * Copyright 2017 <copyright holder> <email>
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mtoc_simple.h"

MTOC_simple::MTOC_simple()
{
    Vector3d position_2( 0.0 , 0.0 , 0.0 );
    this->position = position_2;
    this->friction = 1e-7;
}


MTOC_simple::MTOC_simple( Vector3d position_tmp )
{
    this->position = position_tmp;
    double d_v = 0.5e-6;
    double K = 1.0;
    double effectiveF = sim_of_Cell::multiply_friction_constant * 3.0 * sim_of_Cell::PI * sim_of_Cell::viscosity * d_v * K;
//	this->effectiveFriction = effectiveF; 
    //cout<<"effectiveF = "<<effectiveF<<endl;
    this->friction = 1e-7;
    
    
    
}


MTOC_simple::MTOC_simple(const MTOC_simple& other)
{
    this->position = other.position;
    this->friction = other.friction;
}

Vector3d MTOC_simple::get_position()
{
    return this->position;
}

void MTOC_simple::set_position( Vector3d position_tmp )
{
    this->position = position_tmp;
}

double MTOC_simple::get_friction()
{
    return this->friction;
}

void MTOC_simple::set_friction( double friction_tmp )
{
    this->friction = friction_tmp;
}

void MTOC_simple::oneStepMidStepAlgorithm_1_half( Vector3d force )
{
    double timeHalfStep =  sim_of_Cell::time_Step / 2.0;
    MatrixXd V_0 = force / this->friction; 
	MatrixXd R_half = this->position + V_0 * timeHalfStep;
	this->set_position( R_half );
}


void MTOC_simple::oneStepMidStepAlgorithm_2_half( Vector3d force , Vector3d position_tmp )
{
    MatrixXd V_0 = ( force ) / this->friction; 
    MatrixXd R_final =  position_tmp + V_0 * sim_of_Cell::time_Step;
    this->set_position( R_final );
}







MTOC_simple::~MTOC_simple()
{

}

MTOC_simple& MTOC_simple::operator=(const MTOC_simple& other)
{
    this->position = other.position;
    this->friction = other.friction;  
    return *this;
}
/*
bool MTOC_simple::operator==(const MTOC_simple& other) const
{

}
*/
