/*  This file is part of Victor.

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
/*!
 *  \author    Luca Demo
 *  \date      2015
 */
#include <MobiProtein.h>
#include <PdbSaver.h>

using namespace Victor;
using namespace Victor::Mobi;
using namespace Victor::Biopool;

/**
 * Load Pdb data into Protein object. The specified models are loaded.
 * There is no direct correlation between models id in pdb file and model
 * index inside protein object. Index inside protein depends on the order
 * in which the models are loaded. *
 * @param pl (PdbLoader&) reference to the PdbLoader to load from
 * @param models (vector<unsigned int>) models to load (accordind to pdb file model names),
 */
vector<unsigned int> MobiProtein::load(PdbLoader& pl, vector<unsigned int> models){
	if (verbose)
		cout << "Loading protein models..." <<endl;

	if (!verbose)
		pl.setNoVerbose();
	//Load models
	for(unsigned int i = 0; i < models.size(); i++){
		if (verbose)
			cout << "\t>>>model#" << models[i] << endl;
		pl.setModel(models[i]);
		pl.checkModel();

		this->Protein::load(pl);
		this->modelsID.push_back(models[i]);
	}
	return models;
}

/**
 * Load Pdb data into Protein object. All models are loaded.
 * There is no direct correlation between models id in pdb file and model
 * index inside protein object. Index inside protein depends on the order
 * in which the models are loaded. *
 * @param pl (PdbLoader&) reference to the PdbLoader to load from
 */
vector<unsigned int> MobiProtein::load(PdbLoader& pl){
	vector<unsigned int> models;
	//Load Models
	for(unsigned int i = 1; i <= pl.getMaxModels(); i++)
		models.push_back(i);
	return load(pl,models);
}

/**
 * Get the Spacer object representing the model of this Protein
 * @param model (unsigned int) the model to get
 * @return (Spacer&) the model
 */
Spacer* MobiProtein::getModel(unsigned int model){
	return this->getSpacer(model);
}
