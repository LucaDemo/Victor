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
 * @file PreteinModel.h
 * @author Luca Demo
 * @date May 2015
 * @version 0.1
 */
#ifndef MOBI_SOURCES_PROTEINMODEL_
#define MOBI_SOURCES_PROTEINMODEL_

#include<Protein.h>
#include<PdbLoader.h>
#include<AtomCode.h>

using namespace Victor::Biopool;

namespace Victor { namespace Mobi {
	/**
	 * @brief Extends Protein class with functionalities related to manipulation of NMR model and models comparations.
	 * Scale Distance metric is provided.
	 */
	class ProteinModel : public Protein{
	public:
		ProteinModel() : Protein(), verbose(false){modelsID = vector<unsigned int>(0);};
		ProteinModel(const Protein& _orig) : Protein(_orig), verbose(false){};

		ProteinModel* clone(){
			ProteinModel* tmp = new ProteinModel();
			tmp->copy(*this);
			return tmp;
		}

		void copy(const ProteinModel& orig) {
		    Protein::copy(orig);
		}


		vector<unsigned int> load(PdbLoader& pl, vector<unsigned int> models);
		vector<unsigned int> load(PdbLoader& pl);

		Spacer* getModel(unsigned int _model);
		void SD(vector<double>& sd, ProteinModel& ref, AtomCode atom, double d0 = 4);
		void SD(vector<double>& sd, Spacer* ref, AtomCode atom, double d0 = 4);

		void setVerbose(bool v){
			verbose = v;
		}

		const vector<unsigned int>& getModelsID(){ return modelsID;}

	private:
		bool verbose;
		vector<unsigned int> modelsID;
	};
}}


#endif /* MOBI_SOURCES_PROTEINMODEL_ */
