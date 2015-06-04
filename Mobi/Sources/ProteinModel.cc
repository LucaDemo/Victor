/* 	This file is part of Victor

    Victor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Victor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file ProteinModel.cc
 * @author Luca Demo
 * @date May 2015
 * @version 0.1
 */
#include <ProteinModel.h>
#include <PdbSaver.h>

using namespace Victor;
using namespace Victor::Mobi;
using namespace Victor::Biopool;

/**
 * Load Pdb data into Protein object with the specified model
 * @param pl (PdbLoader&) reference to the PdbLoader to load from
 * @param chain (char) name of the chain to load, * = first chain (default)
 * @param model (unsigned int) model to load, 0 = all models (default)
 */
void ProteinModel::load(PdbLoader& pl, char chain, unsigned int model){
	if (chain != '*')
		pl.setChain(chain);
	else
		pl.setChain(pl.getAllChains()[0]);
	pl.checkAndSetChain();
	if (model > 0){	//load a single model (the specified one)
		pl.setModel(model);
		pl.checkModel();
		this->Protein::load(pl);
	}else			//load all models
		for (unsigned int i = 0; i < pl.getMaxModels(); i++){
			pl.setModel(i);
			pl.checkModel();
			this->Protein::load(pl);
		}
}

/**
 * Get the Spacer object representing the model of this Protein
 * @param model (unsigned int) the model to get
 * @return (Spacer&) the model
 */
Spacer& ProteinModel::getModel(unsigned int model){
	return *(this->getSpacer(model));
}

/**
 * Get Scale Distance between this model and the one in the specified ProteinModel
 * @param ref (ProteinModel&) the model to compare
 * @param atoms (vector<string>&) names of the atoms to consider
 * @param d0 (double) normalization parameter, 4 is default
 */
void ProteinModel::SD(vector<double>& sd, ProteinModel& ref, AtomCode atom, double d0){
	ProteinModel::SD(sd, ref.getModel(0), atom, d0);
}

/**
 * Get Scale Distance between this  model and the one specified
 * @param ref (Spacer&) the model to compare
 * @param atoms (vector<string>&) names of the atoms to consider
 * @param d0 (double) normalization parameter, 4 is default
 */
void ProteinModel::SD(vector<double>& sd, Spacer& ref, AtomCode atom, double d0){
	Spacer& model = this->getModel(0);
	sd.clear();
	if (ref.size() != model.size())
		ERROR("Reference protein spacer has a different number of Amino",exception);

	for (unsigned int i = 0; i < ref.size(); i++){	//foreach Amino
		if (ref.getAmino(i).getType() != model.getAmino(i).getType())
			ERROR("Reference protein Aminos are not compatible",exception);
		Atom mAtom = model.getAmino(i).getAtom((unsigned int)atom);
		Atom refAtom = ref.getAmino(i).getAtom((unsigned int)atom);

		double dist;
		dist = sqrt(pow(mAtom.getCoords().x - refAtom.getCoords().x,2.0)
				+ pow(mAtom.getCoords().y - refAtom.getCoords().y,2.0)
				+ pow(mAtom.getCoords().z - refAtom.getCoords().z,2.0)
		);
		sd.insert(sd.end(), dist);

	}
}


