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
	 * @param _d0 (double) d0 normalization value for scaled distance
	 * @param _atom (AtomCode) atom to use for distance calculations
	 * @param _sd_th (double) Scaled distance mean threshold
	 * @param _sdsd_th (double) Scaled distance deviation threshold
	 * @param _psi_th (double) Psi angles deviation threshold
	 * @param _phi_th (double) Phi angles deviation threshold
	 */
	MobiMethods(double _d0 = DEF_D0, AtomCode _atom = DEF_ATOM, double _sd_th = DEF_SD_TH,
			double _sdsd_th = DEF_SDSD_TH, double _psi_th = DEF_PSI_TH, double _phi_th = DEF_PHI_TH) :
			d0(_d0), atom(_atom), sd_th(_sd_th), sdsd_th(_sdsd_th), psi_th(_psi_th), phi_th(_phi_th), verbose(DEF_VERBOSE){}

	/**
	 * Set verbosity level
	 * @param v (int) verbosity to set
	 */
	void verbosity(int v){
		this->verbose = v;
	}

	/**
	 * Get verbosity level
	 * @return (int) verbosity
	 */
	int verbosity(){
		return this->verbose;
	}

	/**
	 * Calculates the Scaled and "simple" distances given the ProteinModel object containing all the models.
	 * For all possible pair of models, Superimposition is performed and distances are calculated.
	 * @param protein (ProteinModel&) protein object containing all the models
	 * @param tm (TMScoreBin&) TMScore binary binder
	 * @param scaledDist (VectoCollection&) destination of scaled distances
	 * @param dist (VectoCollection&) destination of distances
	 * @param mm (MobiMethods&) parameters
	 */
	static void distances(ProteinModel* protein, TMScoreBinder* tm, VectorCollection& scaledDist, VectorCollection& dist, MobiMethods &mm);


	/**
	 * Calculate the scaled distance given two spacers
	 * @param mod1 (Spacer&) first model
	 * @param mod2 (Spacer&) second model
	 * @param atom (AtomCode) atom to base distance calculation on, default atom if not specified (CA)
	 * @param d0 (double) d0 normalization value, default value if not specified (4)
	 */
	static vector<double> scaledDistance(Spacer* mod1, Spacer* mod2, AtomCode atom = DEF_ATOM, double d0 = DEF_D0);

	/**
	 * Calculate the "simple" distance given two spacers.
	 * @param mod1 (Spacer&) first model
	 * @param mod2 (Spacer&) second model
	 * @param atom (AtomCode) atom to base distance calcutation on, default atom if not specified (CA)
	 */
	static vector<double> distance(Spacer* mod1, Spacer* mod2, AtomCode atom = DEF_ATOM);

	/**
	 * Given a protein containing models, populate a VectorCollection of phi angles
	 * @param protein (ProteinModel*) the protein to elaborate
	 * @param phis (VectorCollection&) output VectorCollection
	 */
	static void phis(ProteinModel* protein, VectorCollection& phis);

	/**
	 * Given a protein containing models, populate a VectorCollection of psi angles
	 * @param protein (ProteinModel*) the protein to elaborate
	 * @param psis (VectorCollection&) output VectorCollection
	 */
	static void psis(ProteinModel* protein, VectorCollection& psis);

	/**
	 * Given a protein containing models, build DSSP mobility estimation.\n The returning
	 * vector contains values according to Mobi Methods guidelines:\n
	 * 0 = non mobile\n
	 * 1 = mobile\n
	 * 2 = coil
	 * @param protein (ProteinModel*) the protein
	 * @return (vector<int>) a vector containing the DSSP estimation (0,1,2 values)
	 */
	static vector<int> DSSP(ProteinModel* protein);


	static vector<int> mobiMobility(ProteinModel* protein, TMScoreBinder* tm, MobiMethods& settings);

	static vector<int> SDFilters(vector<int> sdm, vector<int> const &sdsd,
			vector<int> const &dssp, vector<int> const &phis, vector<int> const &psis, MobiMethods &settings);

	//static ProteinModel* averageModel(ProteinModel* protein, TMScoreBinder* tm);
	double getD0(){
		return this->d0;
	}

	AtomCode getAtom(){
		return this->atom;
	}

	double getSDTh(){
		return this->sd_th;
	}

	double getSDSDTh(){
		return this->sdsd_th;
	}

	double getPhiTh(){
		return this->phi_th;
	}

	double getPsiTh(){
		return this->psi_th;
	}

private:
	double d0;
	double sdsd_th;
	double sd_th;
	double psi_th;
	double phi_th;
	AtomCode atom;
	int verbose;
};

}}	//Namespaces


#endif /* MOBI_SOURCES_MOBIMETHODS_H_ */
