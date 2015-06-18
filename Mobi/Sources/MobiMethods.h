/*
 * MobiMethods.h
 *
 *  Created on: 18/giu/2015
 *      Author: luca
 */

#ifndef MOBI_SOURCES_MOBIMETHODS_H_
#define MOBI_SOURCES_MOBIMETHODS_H_

#include <Atom.h>
#include <Spacer.h>

using namespace Victor::Biopool;
//using namespace Victor::Mobi;

namespace Victor{ namespace Mobi{

extern const double DEF_D0;
extern const AtomCode DEF_ATOM;


class MobiMethods{
public:

	/**
	 * Calulates the average Scaled distance of this results
	 */
	static void avgScaledDistance();


	/**
	 * Calculate the scaled distance given two spacers (previously superimposed with TMScore)
	 * @param sd (vector<double>&) reference to sd output vector
	 * @param mod (Spacer&) model spacer
	 * @param ref (Spacer&) reference spacer
	 * @param atom (AtomCode) atom to base distance calculation on, default atom if not specified
	 * @param d0 (double) d0 normalization value, default value if not specified
	 */
	static void scaledDistance(vector<double>& sd, Spacer& mod, Spacer& ref, AtomCode atom, double d0);
};

}}


#endif /* MOBI_SOURCES_MOBIMETHODS_H_ */
