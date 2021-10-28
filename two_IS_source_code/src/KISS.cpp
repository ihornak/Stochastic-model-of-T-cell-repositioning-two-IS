/*
 * KISS.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: hornak
 */
#include "KISS.h"

/*
namespace KISS
{
*/



void seed_generator_with_time()
{

	//srand(time(NULL)); // use current time as seed for random generator
        unsigned int random_variable_1 = genrand64_int64();
	unsigned int random_variable_2 = genrand64_int64();
        unsigned int random_variable_3 = genrand64_int64();
	unsigned int random_variable_4 = genrand64_int64();

	x = random_variable_1;
	c = random_variable_2;
	y = random_variable_3;
	z = random_variable_4;

	cout<<"x = "<<x<<endl;
	cout<<"c = "<<c<<endl;
	cout<<"y = "<<y<<endl;
	cout<<"z = "<<z<<endl;

}



double KISS_function()
{
	return  ( double )KISS;
}



double randx()
{
	return ( double )KISS / ( double ) ( ULLONG_MAX  + 1.0);
	//return ( double )genrand64_int64() / ( double ) ( ULLONG_MAX  + 1.0);
}

double rand_x( double min ,  double max )
{
	return randx() * ( max - min ) + min;
}

double normalRandom()
{
  double u1=randx();
  double u2=randx();
  return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}

double normal_Random( double mu, double sigma)
{
  double ran = mu + sigma * normalRandom();
  return ran;
}

int RAND_INT( int min , int max )
{
	//this will return integer from min to max - 1 
	double lower = ( double ) min;
	double upper = ( double ) max - 1e-12; // boudary
	double tmp = rand_x( lower ,  upper );
//	cout<<tmp<<endl;
	int final = ( unsigned int )  tmp;
	return final;
}



double triangular_distribution()
{
	double R = IS_Cortical_Sl_parameter::radius;
	double x = randx();
	
	double result = sqrt( x ) * R;
	//cout<<"result = "<<result<<endl;	
	return result;
}




//}



















