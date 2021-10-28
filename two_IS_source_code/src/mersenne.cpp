#include "mersenne.h"


namespace mersenne_twisters_64
{
	//The number of random-generated-numbers that is discarded due to the warm-up
	unsigned int discard_number = 10000;
	//double ext_number_of_generators = 1;
	std::vector< std::mt19937_64 > vector_of_mersenne_twisters;
	std::vector<int> seeds;
	void initialize_generators( unsigned int number_of_threads )
	{
	        //Seed based on time	
		unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
		// Creates generator that will seed other generators used in the simulation
		std::mt19937_64 seeding_generator( seed1 );

		for( unsigned int index = 0 ; index < number_of_threads ; index ++ )
		{
			unsigned int seed_tmp = seeding_generator();
			seeds.push_back( seed_tmp );
			std::mt19937_64 generator_for_array( seed_tmp ); //discard
			//First discard_number of numbers are thrown away due to the warm-up			
			generator_for_array.discard( discard_number );
			//Prepared generator is inserted to the vector of generators			
			vector_of_mersenne_twisters.push_back( generator_for_array );
		}
		//The seeds of generators are stored		
		char name_of_MTOC_file [55];
	  	sprintf ( name_of_MTOC_file , "picturesVideos/numerical_results/seeds.txt" );
	  	FILE *out_MTOC;
	  	out_MTOC = fopen( name_of_MTOC_file , "w");
       
	  	for( unsigned int index = 0 ; index < number_of_threads ; index ++ )
	  	{
			  fprintf( out_MTOC , "%u\n",  seeds[index] );

	  	}
	 	fclose( out_MTOC );
	}


	//This function initializes the generator with the stable seed - it is unused
	void initialize_generators_stable_seed( unsigned int number_of_threads )
	{
		unsigned int seed1 = 4568781;
		std::mt19937_64 seeding_generator(seed1);
		for( unsigned int index = 0 ; index < number_of_threads ; index ++ )
		{
			unsigned int seed_tmp = seeding_generator();
			std::mt19937_64 generator_for_array( seed_tmp ); //discard
			generator_for_array.discard( discard_number );
			vector_of_mersenne_twisters.push_back( generator_for_array );
		}



	}
	
        //This is the generation of the number from the triangular distribution between 0 and the upper_boundary
	double triangular_distribution( double upper_boundary )
	{

		std::uniform_real_distribution<> distribution{ 0 , 1 };
		unsigned int number_of_generator = omp_get_thread_num();
		double R = upper_boundary;
		double x = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
			    
		double result = sqrt( x ) * R;
		return result;
	}







}





void init_genrand64(unsigned long long seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++)
        mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

void init_by_array64(unsigned long long init_key[],
		     unsigned long long key_length)
{
    unsigned long long i, j, k;
    init_genrand64(19650218ULL);
    i=1; j=0;
    k = (NN>key_length ? NN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845ULL))
          + init_key[j] + j;
        i++; j++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 2862933555777941757ULL))
          - i;
        i++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
    }

    mt[0] = 1ULL << 63;
}


unsigned long long genrand64_int64(void)
{
    int i;
    unsigned long long x;
    static unsigned long long mag01[2]={0ULL, MATRIX_A};

    if (mti >= NN) {
        if (mti == NN+1)
            init_genrand64(5489ULL);

        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

        mti = 0;
    }

    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
}



double randx_mersenne()
{
	return ( double )genrand64_int64() / ( double ) ( ULLONG_MAX  + 1.0);
}


double rand_x_mersenne( double min ,  double max )
{
	return randx_mersenne() * ( max - min ) + min;
}


double normalRandom_mersenne()
{
  double u1=randx_mersenne();
  double u2=randx_mersenne();
  return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}

int RAND_INT_mersenne( int min , int max )
{

	double lower = ( double ) min;
	double upper = ( double ) max;
	double tmp = rand_x_mersenne( lower ,  upper );

	int final = ( unsigned int )  tmp;
	return final;
}
