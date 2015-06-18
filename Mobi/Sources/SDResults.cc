/*
 * MobiMethods.cc
 *
 *  Created on: 09/giu/2015
 *      Author: luca
 */


#include <SDResults.h>
using namespace Victor;
using namespace Victor::Biopool;
using namespace Victor::Mobi;

double const DEF_D0 = 4;
AtomCode const DEF_ATOM = CA;

/**
 * Get Scale Distance between two models
 * @param sd (vector<double>&) output vector
 * @param mod (Spacer&) model structure
 * @param ref (Spacer&) reference structure
 * @param atoms (vector<string>&) names of the atoms to consider, default value is used if not specified
 * @param d0 (double) normalization parameter, default value is used if not specified
 */
void Victor::Mobi::scaledDistance(vector<double>& sd, Spacer& mod, Spacer& ref, AtomCode atom, double d0){
	sd.clear();
	if (ref.size() != mod.size())
		ERROR("Reference protein spacer has a different number of Amino",exception);
	for (unsigned int i = 0; i < ref.size(); i++){	//foreach Amino
		if (ref.getAmino(i).getType() != mod.getAmino(i).getType())
			ERROR("Reference protein Aminos are not compatible",exception);
		Atom mAtom = mod.getAmino(i).getAtom((unsigned int)atom);
		Atom refAtom = ref.getAmino(i).getAtom((unsigned int)atom);

		double dist;
		dist = sqrt(pow(mAtom.getCoords().x - refAtom.getCoords().x,2.0)
				+ pow(mAtom.getCoords().y - refAtom.getCoords().y,2.0)
				+ pow(mAtom.getCoords().z - refAtom.getCoords().z,2.0)
		);
		sd.insert(sd.end(), dist);
	}
}

/**
 * Calculates the mean Scaled Distance by calculating the mean of all
 * the distance between pair of aligned models.
 * @param mean (vector<double> *) pointer to output vector
 */
void SDResults::meanSD(vector<double> *mean){
	mean = new vector<double>(this->modelsSize());
	std::map<int,std::vector<double> >::const_iterator it;
	for (it = this->results->begin(); it != this->results->end(); ++it)
		for (unsigned int a = 0; a < this->modelsSize(); a++)
			(*mean)[a] += it->second[a];
	for (unsigned int a = 0; a < this->modelsSize(); a++)
		(*mean)[a] /= this->size();
}
