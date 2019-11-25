/*
 * math_functions.h
 *
 *  Created on: Oct 7, 2019
 *      Author,Copyright: Rene Staritzbichler
 */

#ifndef MATH_FUNCTIONS_H_
#define MATH_FUNCTIONS_H_



float SquaredDistance( const std::vector< float> & P1, const std::vector< float > & P2)
{
	float d = 0.0;

	for( std::vector< float >::const_iterator i1 = P1.begin(), i2 = P2.begin(); i1 != P1.end(); ++i1, ++i2)
    {
		d += pow( *i1 - *i2 , 2);
    }

	return d;
}



float RMSD( const std::vector< std::vector< float > > & P1, const std::vector< std::vector< float > > & P2)
{
	if( P1.size() != P2.size() || P1.size() == 0)
    {
		std::cout << "ERROR: mismatch in " << __FUNCTION__ << " " << P1.size() << " " << P2.size() << std::endl;
		exit(1);
    }

	float rmsd = 0.0;
	for( std::vector< std::vector< float > >::const_iterator i1 = P1.begin(), i2 = P2.begin(); i1 != P1.end(); ++i1, ++i2)
    {
		rmsd += SquaredDistance( *i1, *i2);
    }

	return sqrt( rmsd / (float) P1.size() );
}



#endif /* MATH_FUNCTIONS_H_ */
