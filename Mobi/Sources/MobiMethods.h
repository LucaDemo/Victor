/*
 * MobiMethods.h
 *
 *  Created on: 18/giu/2015
 *      Author: luca
 */

#ifndef MOBI_SOURCES_MOBIMETHODS_H_
#define MOBI_SOURCES_MOBIMETHODS_H_

#include <Atom.h>
#include <AtomCode.h>
#include <Spacer.h>
#include <VectorCollection.h>
#include <TMScoreBin.h>

using namespace Victor::Biopool;
//using namespace Victor::Mobi;

namespace Victor{ namespace Mobi{

extern const double DEF_D0;
extern const AtomCode DEF_ATOM;

extern const int MOBI_VERBOSE_LEVEL;


class MobiMethods{
public:

	/**
	 * Constructor with parameters
	 * @param _d0 (double) d0 normalization value for scaled distance
	 * @param _atom (AtomCode) atom to use for distance calculations
	 */
	MobiMethods(double _d0 = DEF_D0, AtomCode _atom = DEF_ATOM) : d0(_d0), atom(_atom){}

	/**
	 * Calculates the Scaled distances given the ProteinModel object containing all the models.
	 * For all possible pair of models, Superimposition is performed and scaled distance is recorded.
	 * @param protein (ProteinModel&) protein object containing all the models
	 * @param tm (TMScoreBin&) TMScore binary binder
	 * @return (VectoCollection*) reference to an object containing all the distances
	 */

	static VectorCollection* scaledDistances(ProteinModel& protein, TMScoreBin& tm);


	/**
	 * Calculate the scaled distance given two spacers (previously superimposed with TMScore)
	 * @param mod (Spacer&) model spacer
	 * @param ref (Spacer&) reference spacer
	 * @param atom (AtomCode) atom to base distance calculation on, default atom if not specified
	 * @param d0 (double) d0 normalization value, default value if not specified
	 */
	static vector<double> scaledDistance(Spacer& mod1, Spacer& mod2, AtomCode atom = DEF_ATOM, double d0 = DEF_D0);

	/**
	 * Calculate the scaled distance given two spacers (previously superimposed with TMScore).
	 * Parameters passed as reference to MobiMethods object
	 * @param mod (Spacer&) model spacer
	 * @param ref (Spacer&) reference spacer
	 * @param pars (MobiMethod&) reference to MobiMethod instance with parameters
	 */
	static vector<double> scaledDistance(Spacer& mod1, Spacer& mod2, MobiMethods& pars){
		return scaledDistance(mod1, mod2, pars.getAtom(), pars.getD0());
	}

	double getD0(){
		return this->d0;
	}

	AtomCode getAtom(){
		return this->atom;
	}
private:
	double d0;
	AtomCode atom;
};

}}


#endif /* MOBI_SOURCES_MOBIMETHODS_H_ */
