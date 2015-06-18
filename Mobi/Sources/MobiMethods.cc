/*
 * MobiMethods.cc
 *
 *  Created on: 09/giu/2015
 *      Author: luca
 */

#include <Atom.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <MobiMethods.h>

using namespace Victor::Biopool;
using namespace Victor::Mobi;


/**
 * Default d0 value: 4 as specified in "Mobi Methods" documentation
 */
double const DEF_D0 = 4;
/**
 * Default atom: Alpha Carbon
 */
AtomCode const DEF_ATOM = CA;

void MobiMethods::avgScaledDistance(){}

void MobiMethods::scaledDistance(vector<double>& sd, Spacer& mod, Spacer& ref, AtomCode atom, double d0){
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
