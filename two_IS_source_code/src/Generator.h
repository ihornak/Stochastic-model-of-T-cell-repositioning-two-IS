/*
 * generator.h
 *
 *  Created on: Jul 24, 2015
 *      Author: hornak
 */

#ifndef GENERATOR_H_
#define GENERATOR_H_
#include <random>
#include <iostream>


//This is very basic random number generator
//Call only if the generators for parallel processing are overkill

class Generator {
public:


	Generator( unsigned int seed ); 
	double get_probability();
	


	virtual ~Generator();


private:

	std::mt19937_64 rd;
	std::uniform_real_distribution<> dist_uniform;
    
};

#endif /* MICROTUBULE_H_ */
