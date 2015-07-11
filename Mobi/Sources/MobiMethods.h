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
extern const double DEF_PHI_TH;
extern const double DEF_PSI_TH;
extern const double DEF_SD_TH;
extern const double DEF_SDSD_TH;

extern const int DEF_VERBOSE;

/**
 * This class provides static methods to calculate mobility with Mobi methods.\n
 * Can be instantiated to provide settings to the static methods. Default constructor provides
 * default mobi parameters.
 */
class MobiMethods{
public:

	/**
	 * Constructor with parameters.
	 * @param _d0  d0 normalization value for scaled distance
	 * @param _atom  atom to use for distance calculations
	 * @param _sd_th  Scaled distance mean threshold
	 * @param _sdsd_th  Scaled distance deviation threshold
	 * @param _psi_th  Psi angles deviation threshold
	 * @param _phi_th  Phi angles deviation threshold
	 */
	MobiMethods(double _d0 = DEF_D0, AtomCode _atom = DEF_ATOM, double _sd_th = DEF_SD_TH,
			double _sdsd_th = DEF_SDSD_TH, double _psi_th = DEF_PSI_TH, double _phi_th = DEF_PHI_TH) :
			 d0(_d0), atom(_atom), sd_th(_sd_th), sdsd_th(_sdsd_th), psi_th(_psi_th), phi_th(_phi_th), verbose(DEF_VERBOSE){}

	/**
	 * Set verbosity level
	 * @param v  verbosity to set
	 */
	void verbosity(int v){
		this->verbose = v;
	}

	/**
	 * Get verbosity level
	 * @return verbosity
	 */
	int verbosity(){
		return this->verbose;
	}

	/**
	 * Calculates the Scaled and "simple" distances given the ProteinModel object containing all the models.
	 * For all possible pair of models, Superimposition is performed and distances are calculated.
	 * @param protein  pointer to protein object containing all the models
	 * @param tm pointer to TMScore binary binder
	 * @param scaledDist destination of scaled distances
	 * @param dist destination of distances
	 * @param mm reference to MobiMethod object for mobi parameters
	 */
	static void distances(ProteinModel* protein, TMScoreBinder* tm, VectorCollection& scaledDist, VectorCollection& dist, MobiMethods &mm);


	/**
	 * Calculate the scaled distance given two models
	 * @param mod1 the first model (Spacer pointer)
	 * @param mod2 the second model (Spacer pointer)
	 * @param atom atom to base distance calculation on, default atom if not specified (CA)
	 * @param d0  d0 normalization value, default value if not specified (4)
	 */
	static vector<double> scaledDistance(Spacer* mod1, Spacer* mod2, AtomCode atom = DEF_ATOM, double d0 = DEF_D0);

	/**
	 * Calculate the "simple" distance given two models.
	 * @param mod1 the first model (Spacer pointer)
	 * @param mod2 the second model (Spacer pointer)
	 * @param atom (AtomCode) atom to base distance calcutation on, default atom if not specified (CA)
	 */
	static vector<double> distance(Spacer* mod1, Spacer* mod2, AtomCode atom = DEF_ATOM);

	/**
	 * Given a protein containing models, populate a VectorCollection of phi angles
	 * @param protein pointer to the protein to elaborate
	 * @param phis output VectorCollection
	 */
	static void phis(ProteinModel* protein, VectorCollection& phis);

	/**
	 * Given a protein containing models, populate a VectorCollection of psi angles
	 * @param protein pointer to the protein to elaborate
	 * @param psis output VectorCollection
	 */
	static void psis(ProteinModel* protein, VectorCollection& psis);

	/**
	 * Given a protein containing models, build DSSP mobility estimation.\n The returning
	 * vector contains values according to Mobi Methods guidelines:\n
	 * 0 = non mobile\n
	 * 1 = mobile\n
	 * 2 = coil
	 * @param protein reference to the protein
	 * @return a vector<int> containing the secondary mobility estimation (0,1,2 values)
	 */
	static vector<int> secondaryMobi(ProteinModel* protein);

	/**
	 * Performs mobility calculation
	 * @param protein pointer to the protein to process
	 * @param tm pointer to TM imposer
	 * @param settings reference to MobiMethods object for Mobi parameters
	 * @param output reference to output destination (fasta), default is cout
	 * @return a vector<int> of mobility values (0,1)
	 */
	vector<int> mobiMobility(ProteinModel* protein, TMScoreBinder* tm);

	/**
	 * Applies Mobi filters to Scaled Distance mobility track to obtain final mobility values
	 * @param sdm vector<int> of scaled distance mobility track
	 * @param sdsd vector<int> scaled distance std dev mobility track
	 * @param sec vector<int> secondary structure mobility track
	 * @param phis vector<int> phi angles mobility track
	 * @param psis vector<int> psi angles mobility track
	 * @param settings reference to MobiMethods object for mobi parameters
	 * @return vector<int> of filtered mobility ("all" mobility)
	 */
	static vector<int> SDFilters(vector<int> const &sdm, vector<int> const &sdsd,
			vector<int> const &sec, vector<int> const &phis, vector<int> const &psis, MobiMethods &settings);

	/**
	 * Applies Mobi filters to Scaled Distance mobility track to obtain final mobility values
	 * @param sdm vector<int> of scaled distance mobility track
	 * @param sdsd vector<int> scaled distance std dev mobility track
	 * @param sec vector<int> secondary structure mobility track
	 * @param phis vector<int> phi angles mobility track
	 * @param psis vector<int> psi angles mobility track
	 * @param settings reference to MobiMethods object for mobi parameters
	 * @return vector<int> of filtered mobility ("all" mobility)
	 */
	static vector<int> SDFilters(MobiMethods &settings){
		return SDFilters(settings.getSDMeanMobility(), settings.getSDDevsMobility(), settings.getSecMobility(), settings.getPhiMobility(),
				settings.getPsiMobility(), settings);
	}



	//static ProteinModel* averageModel(ProteinModel* protein, TMScoreBinder* tm);
	double getD0(){	return this->d0;}

	AtomCode getAtom(){return this->atom;}

	double getSDTh(){return this->sd_th;}

	double getSDSDTh(){return this->sdsd_th;}

	double getPhiTh(){return this->phi_th;}

	double getPsiTh(){return this->psi_th;}


	const vector<double> getSDMeans(){return SDMeans;}

	const vector<int> getSDMeanMobility(){return SDMeanMobility;}

	const vector<double> getSDDevs(){return SDDevs;}

	const vector<int> getSDDevsMobility(){return SDDevsMobility;}

	const vector<double> getPsisAngles(){return psisAngles;}

	const vector<int> getPsiMobility(){return PsiMobility;}

	const vector<double> getPhisAngles(){return phisAngles;}

	const vector<int> getPhiMobility(){return PhiMobility;}

	const vector<int> getSecMobility(){return secMobility;}

	const vector<int> getMobiMobility(){return mobiMob;}

	const vector<char> getSequence(){return sequence;}

	VectorCollection& getDistances(){return dist;}

	VectorCollection& getScaledDistances(){return scaledDist;}


private:
	double d0;
	double sdsd_th;
	double sd_th;
	double psi_th;
	double phi_th;
	AtomCode atom;
	int verbose;

	vector<char> sequence;
	vector<double> SDMeans;
	vector<int> SDMeanMobility;
	vector<double> SDDevs;
	vector<int> SDDevsMobility;
	vector<double> psisAngles;
	vector<int> PsiMobility;
	vector<double> phisAngles;
	vector<int> PhiMobility;
	vector<int> secMobility;
	vector<int> mobiMob;

	VectorCollection scaledDist;	//Scaled distances
	VectorCollection dist;			//"Simple" distances
};

}}	//Namespaces


#endif /* MOBI_SOURCES_MOBIMETHODS_H_ */
