/*
 * Nucleus_parametres.h
 *
 *  Created on: Jan 25, 2017
 *      Author: hornak
 */

#ifndef NUCLEUS_PARAMETRES_H_
#define NUCLEUS_PARAMETRES_H_


namespace Nucleus_parametres
{
	//Parameters used for the calculation of the force of the nucleus
	const double wall_nucleus_k1 = 6e-6;//6e-6
	const double wall_nucleus_k2 = 2.0;//
	//Axis of the nucleus
	const double A_AXIS = 3.8e-6;
	const double B_AXIS = 3.8e-6;
	//describes the movement of the nucleus on the z-axis - movement was not implemented so far
        const double z_coordinate = 0.0e-6;
}





#endif /* NUCLEUS_PARAMETRES_H_ */
