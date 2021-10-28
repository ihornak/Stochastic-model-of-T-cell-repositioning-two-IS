#include "Generator.h"


Generator::Generator( unsigned int seed )
{

	//std::random_device rd_tmp;  //Will be used to obtain a seed for the random number engine
    	std::mt19937_64 gen( seed ); //Standard mersenne_twister_engine seeded with rd()
	this->rd = gen;
        std::uniform_real_distribution<> dis(0.0, 1.0); 
        this->dist_uniform = dis;

} 


double Generator::get_probability()
{
        return  this->dist_uniform( this->rd );
}







Generator::~Generator()
{


}


