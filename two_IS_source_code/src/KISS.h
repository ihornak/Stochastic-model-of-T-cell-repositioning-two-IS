/*
 * KISS.h
 *
 *  Created on: Jul 30, 2015
 *      Author: hornak
 * This is the implementation of the Marsaglia's KISS generator(1993)
 * It is not used in the simulations
 * I keep it here just to control Mersenne-Twister
 */

#ifndef KISS_H_
#define KISS_H_


#include <stdio.h>
#include <climits>
#include <math.h>
#include <iostream>
#include "mersenne.h"
#include "simulationofCell.h"
#include "IS_Cortical_Sl_parameter.h"

using namespace std;





/*
namespace KISS
{
*/

static unsigned long long
//x = 7613131461312ULL , c = 258795254ULL,
//y = 13217149525ULL , z = 9758180ULL , t;


x = 9781641336351611ULL , c = 816813135002ULL,
y = 9361448128545ULL , z = 981641ULL , t;

//x = 321354168451611ULL , c = 5002ULL,
//y = 936145ULL , z = 171564616641ULL , t;

//x = 681361363131ULL , c = 9689898864ULL,
//y = 21321323232ULL , z = 546456566ULL , t;

//x = 159159159ULL , c = 879841961474ULL,
//y = 333669877515ULL , z = 13548978468ULL , t;


//x = 9815616141151221ULL , c = 784744152555ULL,
//y = 6987154154111000ULL , z = 325654686548ULL , t;

//x = 23987899519771ULL , c = 987888934126ULL,
//y = 58941376589ULL , z = 9875638191640ULL , t;


//x = 41353584467452ULL , c = 12345678982ULL,
//y = 78548215395ULL , z = 6435158451475ULL , t;


//x = 81961516854ULL , c = 4781616168165ULL,
//y = 96461613542ULL , z = 9816648164158ULL , t;


//x = 46125461346ULL , c = 3132132132333ULL,
//y = 56464546466ULL , z = 9654311156065ULL , t;


//x = 91136161356ULL , c = 1311313231333ULL,
//y = 14789312549ULL , z = 4626245466612ULL , t;


//x = 9843146ULL , c = 66666544885ULL,
//y = 98494684645656ULL , z = 32132133113ULL , t;

//x = 98661616ULL , c = 66544885ULL,
//y = 78986846ULL , z = 2133113ULL , t;


//x = 56456545446ULL , c = 9879898888ULL,
//y = 333214785469689ULL , z = 21331149863ULL , t;


//x = 8768415616513136ULL , c = 114781817181ULL,
//y = 6879831331311313ULL , z = 368366838383ULL , t;

//x = 98983133000306ULL , c = 313321351355ULL,
//y = 6546816513303ULL , z = 9688846611335ULL , t;

	 // use current time as seed for random generator

//static unsigned long long x , c , y , z;



/*
static unsigned long long	x = genrand64_int64();
static unsigned long long	c = genrand64_int64();
static unsigned long long	y = genrand64_int64();
static unsigned long long	z = genrand64_int64();
*/
/*
static unsigned long long t;
static unsigned long long	x;
static unsigned long long	c;
static unsigned long long	y;
static unsigned long long	z;
*/


#define MWC (t=(x<<58)+c, c=(x>>6), x+=t, c+=(x<t), x)
#define XSH ( y^=(y<<13), y^=(y>>17), y^=(y<<43) )
#define CNG ( z=6906969069LL*z+1234567 )
#define KISS (MWC+XSH+CNG)




double KISS_function();

double randx();

double rand_x( double min ,  double max );


double normalRandom();

double normal_Random( double mu, double sigma);
int RAND_INT( int min , int max );

void seed_generator_with_time();


double triangular_distribution();

//}


//









#endif /* KISS_H_ */
