


#include <stdio.h>


#include <stdio.h>
#include <climits>
#include <math.h>
#include <iostream>
#include <random>
#include <chrono>
#include <omp.h>


#ifndef MERSENNE_H_
#define MERSENNE_H_
using namespace std;



#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */



namespace mersenne_twisters_64
{
	//this is the std: vector, where generators are stored
	//every thread has one generator
	extern std::vector< std::mt19937_64 > vector_of_mersenne_twisters;
	// Initialize generators with time 	
	void initialize_generators( unsigned int number_of_threads );
	// Initialize generators with a stable seed 	
	void initialize_generators_stable_seed( unsigned int number_of_threads );
	//this is the generation of the number from the triangular distribution between 0 and the upper_boundary	
	double triangular_distribution( double upper_boundary );



}



//Default functions;
static unsigned long long mt[NN];
static int mti=NN+1;
void init_genrand64(unsigned long long seed);
void init_by_array64(unsigned long long init_key[], unsigned long long key_length);
unsigned long long genrand64_int64(void);
double randx_mersenne();
double rand_x_mersenne( double min ,  double max );
double normalRandom_mersenne();
int RAND_INT_mersenne( int min , int max );



#endif
