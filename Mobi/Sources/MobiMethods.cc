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
#include <MobiUtils.h>
#include <TMScoreBin.h>

using namespace Victor::Biopool;
using namespace Victor::Mobi;

double const Victor::Mobi::DEF_D0 = 4;
AtomCode const Victor::Mobi::DEF_ATOM = CA;
double const Victor::Mobi::DEF_PHI_TH = 20;
double const Victor::Mobi::DEF_PSI_TH = 20;
double const Victor::Mobi::DEF_SD_TH = .85;
double const Victor::Mobi::DEF_SDSD_TH = .09;
int const Victor::Mobi::DEF_VERBOSE = 0;

template<typename T>
/**
 * Helper method for verbose output
 */
void printVector(vector<T> const &v, string pre = ""){
	if (pre != "")
		cout << pre << endl;
	for (unsigned int i = 0; i < v.size(); i++)
		cout << v[i];
	cout << endl;
}

void MobiMethods::distances(ProteinModel* protein, TMScoreBinder* tm, VectorCollection& scaledDist, VectorCollection& dist, MobiMethods &mm){
	scaledDist.clear();
	dist.clear();
	ProteinModel* imposed;
	vector<double> tmpDist;
	//Perform superimposition for every model pair
	for (unsigned int i = 0; i < protein->size(); i++)
		for (unsigned int j = i+1; j < protein->size(); j++){
			tm->TMScore(*protein,i,j,&imposed);
			tmpDist = scaledDistance(imposed->getModel(0),protein->getModel(j), mm.getAtom(), mm.getD0());
			scaledDist.addValue((i*1000 + j), tmpDist);
			tmpDist = distance(imposed->getModel(0),protein->getModel(j), mm.getAtom());
			dist.addValue((i*1000 + j), tmpDist);
			delete imposed;
		}
}

vector<double> MobiMethods::distance(Spacer* mod1, Spacer* mod2, AtomCode atom){
	if (mod2->size() != mod1->size())
		ERROR("Reference protein spacer has a different number of Aminos",exception);
	vector<double> dist(mod1->size());
	for (unsigned int i = 0; i < mod2->size(); i++){	//foreach Amino
		if (mod2->getAmino(i).getType() != mod1->getAmino(i).getType())
			ERROR("Aminos order in the two sequences is not the same",exception);
		Atom& mAtom = mod1->getAmino(i).getAtom(atom);
		Atom& refAtom = mod2->getAmino(i).getAtom(atom);

		dist[i] = sqrt((mAtom.getCoords() - refAtom.getCoords()).square());
	}
	return dist;
}


vector<double> MobiMethods::scaledDistance(Spacer* mod1, Spacer* mod2, AtomCode atom, double d0){
	if (mod2->size() != mod1->size())
		ERROR("Reference protein spacer has a different number of Aminos",exception);
	vector<double> sd(mod1->size());
	for (unsigned int i = 0; i < mod2->size(); i++){	//foreach Amino
		if (mod2->getAmino(i).getType() != mod1->getAmino(i).getType())
			ERROR("Aminos order in the two sequences is not the same",exception);
		Atom& mAtom = mod1->getAmino(i).getAtom(atom);
		Atom& refAtom = mod2->getAmino(i).getAtom(atom);

		double dist = sqrt((mAtom.getCoords() - refAtom.getCoords()).square());
		sd[i] = 1/(1 + pow(dist / d0, 2.0));
	}
	return sd;
}


void MobiMethods::phis(ProteinModel* protein, VectorCollection& phis){
	phis.clear();
	vector<double> modelPhis(protein->getModel(0)->sizeAmino());
	for (unsigned int i = 0; i < protein->size(); i++){
		for (unsigned int j = 0; j < protein->getModel(i)->sizeAmino(); j++){
			AminoAcid& aa = protein->getModel(i)->getAmino(j);
			modelPhis[j] = aa.getPhi();
		}
		phis.addValue(i,modelPhis);
	}
}

void MobiMethods::psis(ProteinModel* protein, VectorCollection& psis){
	psis.clear();
	vector<double> modelPsis(protein->getModel(0)->sizeAmino());
	for (unsigned int i = 0; i < protein->size(); i++){
		for (unsigned int j = 0; j < protein->getModel(i)->sizeAmino(); j++){
			AminoAcid& aa = protein->getModel(i)->getAmino(j);
			modelPsis[j] = aa.getPsi();
		}
		psis.addValue(i,modelPsis);
	}
}

vector<int> MobiMethods::secondaryMobi(ProteinModel* protein){
	//Check for errors in protein
	unsigned int const protLen = protein->size();
	if (protLen < 1)
		ERROR("Empty protein object, no models",exception);
	unsigned int const modLen = protein->getModel(0)->sizeAmino();
	if (modLen < 1)
		ERROR("Empty spacer or incorrect size",exception);
	vector<int> dssp(modLen,0);

	//Build DSSP matrix for this protein
	//1st index = Rows = models; 2nd index = Columns = Residues
	char protDSSP[protLen][modLen];
	for (unsigned int i = 0; i < protLen; i++){	//foreach i model
		//get DSSP for i-model
		vector<set<char> > modelDSSP = protein->getModel(i)->getDSSP();
		if (modelDSSP.size() != modLen)
			ERROR("DSSP length differs from expected",exception);
		for (unsigned int j = 0; j < modLen; j++)	//foreach j amino in i-model
			if (modelDSSP[j].size() > 0)
				protDSSP[i][j] = (*(modelDSSP[j].begin()));
			else
				protDSSP[i][j] = 'C';
	}
	//Build DSSP mobility
	for (unsigned int j = 0; j < modLen; j++){	//foreach residue
		//assign Coil (2) value
		for (unsigned int i = 0; i < protLen; i++) //foreach model
			if (protDSSP[i][j] == 'C' || protDSSP[i][j] == 'S')
				dssp[j] = 2;
			else{
				dssp[j] = 0;
				break;
			}
		//assign Mobile (1) value
		if (dssp[j] != 2)	//not coil
			for (unsigned int i = 0; i < protLen; i++) //foreach model
				if (protDSSP[i][j] != protDSSP[0][j]){
					dssp[j] = 1;
					break;
				}
	}
	return dssp;
}



vector<int> MobiMethods::mobiMobility(ProteinModel* protein, TMScoreBinder* tm){
	//Scaled Distances
    distances(protein,tm,scaledDist,dist,*this);

    //Set sequence
    for (unsigned int i = 0; i < protein->getModel(0)->sizeAmino(); i++)
    	sequence.push_back(protein->getModel(0)->getAmino(i).getType1L());

    //Scaled Distances Mean
    SDMeans = scaledDist.mean();
    SDMeanMobility= vector<int>(SDMeans.size());
    for (unsigned int i = 0; i < SDMeans.size(); i++)
    	SDMeanMobility[i] = SDMeans[i] < sd_th ? 1 : 0;
    //Scaled Distances deviations
    SDDevs = scaledDist.stdDev();
    SDDevsMobility = vector<int>(SDDevs.size());
	for (unsigned int i = 0; i < SDDevs.size(); i++)
		SDDevsMobility[i] = SDDevs[i] > sdsd_th ? 1 : 0;

	VectorCollection angles;
	//Psi angles
	psis(protein,angles);	//values are cleared inside getPhis
	psisAngles = angles.stdDev();
	PsiMobility = vector<int>(psisAngles.size());
	for (unsigned int i = 0; i < psisAngles.size(); i++)
		PsiMobility[i] = psisAngles[i] > psi_th ? 1 : 0;
	//Phi angles
	phis(protein,angles);	//values are cleared inside getPhis
	phisAngles = angles.stdDev();
	PhiMobility = vector<int>(phisAngles.size());
	for (unsigned int i = 0; i < phisAngles.size(); i++)
		PhiMobility[i] = phisAngles[i] > phi_th ? 1 : 0;
	//Secondary mobility
	secMobility = secondaryMobi(protein);

	//Calculate Mobi Mobility
	mobiMob = SDFilters(SDMeanMobility, SDDevsMobility, secMobility, PhiMobility, PsiMobility, *this);

	return mobiMob;
}

vector<int> MobiMethods::SDFilters(vector<int> const &sd, vector<int> const &sdsd,
		vector<int> const &dssp, vector<int> const &phis, vector<int> const &psis, MobiMethods &mm){
	vector<int> sdm = sd;
	if (mm.verbosity() > 1)
		printVector(sdm,"Initial SD");
	//DSSP non-mobile => 1->0 in SD
	for (unsigned int i = 0; i < sdm.size(); i++)
		if (dssp[i] == 0)
			sdm[i] = 0;
	if (mm.verbosity() > 1)
		printVector(sdm,"After DSSP");

	//Patterns in SD
	//1011 -> 1111
	for (unsigned int i = 1; i < sdm.size() - 2; i++)
		if (sdm[i] == 0 && sdm[i-1] == 1 && sdm[i+1] == 1 && sdm[i+2] == 1){
			sdm[i] = 1;
			i+=2;
			if (mm.verbosity() > 1)
				printVector(sdm, "1011->1111");
		}
	//1101 -> 1111
	for (unsigned int i = 2; i < sdm.size() - 1; i++)
		if (sdm[i] == 0 && sdm[i-1] == 1 && sdm[i-2] == 1 && sdm[i+1] == 1){
			sdm[i] = 1;
			i+=1;
			if (mm.verbosity() > 1)
				printVector(sdm, "1101->1111");
		}
	//10011 -> 11111
	for (unsigned int i = 1; i < sdm.size() - 3; i++)
		if (sdm[i] == 0 && sdm[i+1] == 0 && sdm[i-1] == 1 && sdm[i+2] == 1 && sdm[i+3] == 1){
			sdm[i] = 1;
			sdm[i+1] = 1;
			i+=3;
			if (mm.verbosity() > 1)
				printVector(sdm,"10011->11111");
		}
	//11001 -> 11111
	for (unsigned int i = 2; i < sdm.size() -2; i++)
		if (sdm[i] == 0 && sdm[i+1] == 0 && sdm[i-1] == 1 && sdm[i-2] == 1 && sdm[i+2] == 1){
			sdm[i] = 1;
			sdm[i+1] = 1;
			i+=2;
			if (mm.verbosity() > 1)
				printVector(sdm,"11001->11111");
		}
	//01010 -> 00000
	for (unsigned int i = 0; i < sdm.size() - 4; i++)
		if (sdm[i] == 0 && sdm[i+1] == 1 && sdm[i+2] == 0 && sdm[i+3] == 1 && sdm[i+4] == 0){
			sdm[i+1] = 0;
			sdm[i+3] = 0;
			i+=3;
			if (mm.verbosity() > 1)
				printVector(sdm,"01010->00000");
		}
	//00100 -> 00000
	for (unsigned int i = 2; i < sdm.size() - 2; i++)
		if (sdm[i] == 1 && sdm[i+1] == 0 && sdm[i+2] == 0 && sdm[i-1] == 0 && sdm[i-2] == 0){
			sdm[i] = 0;
			i+=2;
			if (mm.verbosity() > 1)
				printVector(sdm,"00100->00000");
		}
	//001100 -> 000000
	for (unsigned int i = 2; i < sdm.size() - 3; i++)
		if (sdm[i] == 1 && sdm[i+1] == 1 && sdm[i-1] == 0 && sdm[i-2] == 0 && sdm[i+2] == 0 && sdm[i+3] == 0){
			sdm[i] = 0;
			sdm[i+1] = 0;
			i+=3;
			if (mm.verbosity() > 1)
				printVector(sdm,"001100->000000");
		}

	//Other patterns
	//110->111
	for (unsigned int i = 0; i < sdm.size() - 2; i++)
		if (sdm[i] == 1 && sdm[i+1] == 1 && sdm[i+2] == 0)
			if (phis[i+2] == 1 && psis[i+2] == 1 && sdsd[i+2] == 1 && psis[i+1] == 1){
				sdm[i+2] = 1;
				if (mm.verbosity() > 1)
					printVector(sdm, "110->111");
			}
	//011->111
	for (unsigned int i = 0; i < sdm.size() - 2; i++)
		if (sdm[i] == 0 && sdm[i+1] == 1 && sdm[i+2] == 1)
			if (phis[i] == 1 && psis[i] == 1 && sdsd[i] == 1 && phis[i+1] == 1){
				sdm[i] = 1;
				if (mm.verbosity() > 1)
					printVector(sdm, "011->111");
			}

	return sdm;
}
