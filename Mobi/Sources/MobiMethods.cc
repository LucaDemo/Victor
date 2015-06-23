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
#include <TMScoreBin.h>

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

VectorCollection* MobiMethods::scaledDistances(ProteinModel& protein, TMScoreBin& tm){
	VectorCollection* distances = new VectorCollection();
	ProteinModel* imposed;
	vector<double> distance;
	//Perform superimposition for every model pair
	for (unsigned int i = 0; i < protein.size(); i++)
		for (unsigned int j = i+1; j < protein.size(); j++){
			cout << "Distance between models: " << i << " and " << j << endl;
//			cout << n << "\tTMI " << i << "-" << j << "...";
			tm.TMImpose(protein,i,j,&imposed);
			distance = MobiMethods::scaledDistance(imposed->getModel(0),protein.getModel(j),CA,4);
//			cout << "OK!"<< endl;
			distances->addValue((i*1000 + j), distance);
			delete imposed;
		}
	//Copy average scaled distance and destroy collection
	return distances;
}


vector<double> MobiMethods::scaledDistance(Spacer& mod1, Spacer& mod2, AtomCode atom, double d0){
	if (mod2.size() != mod1.size())
		ERROR("Reference protein spacer has a different number of Aminos",exception);
	vector<double> sd(mod1.size());
	for (unsigned int i = 0; i < mod2.size(); i++){	//foreach Amino
		if (mod2.getAmino(i).getType() != mod1.getAmino(i).getType())
			ERROR("Aminos order in the two sequences is not the same",exception);
		Atom& mAtom = mod1.getAmino(i).getAtom(atom);
		Atom& refAtom = mod2.getAmino(i).getAtom(atom);

		double dist = sqrt((mAtom.getCoords() - refAtom.getCoords()).square());
		sd[i] = 1/(1 + pow(dist / d0, 2.0));
	}
	return sd;
}
