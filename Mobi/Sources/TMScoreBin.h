/*
 * TMScoreBin.h
 *
 *  Created on: 01/giu/2015
 *      Author: luca
 */

#ifndef MOBI_SOURCES_TMSCOREBIN_H_
#define MOBI_SOURCES_TMSCOREBIN_H_


#include <iostream>
#include <string>
#include <ProteinModel.h>
#include <TMScoreBinder.h>

using namespace Victor::Biopool;
using namespace Victor::Mobi;
using namespace std;


extern const std::string TMTMP_IN1;
extern const std::string TMTMP_IN2;
extern const std::string TMTMP_OUT;


namespace Victor { namespace Mobi {
	/**
	 * @brief TMScore functionalities through external binary.
	 */
	class TMScoreBin: public TMScoreBinder{
	public:

		/**
		 * @brief Setup.
		 * @param _binary (string) full path to binary TMScore file, must have execution permission
		 * @param _tmp (string) full path to temp dir, must have write permission
		 * @param _verbose (bool) verbosity
		 */
		TMScoreBin(std::string _binary = "TMScore", std::string _tmp = ".", bool _verbose = false) :
			binary(_binary),
			tmp(_tmp.substr(_tmp.length()-1,1) == "/" ? _tmp : _tmp + "/"),
			verbose(_verbose)
		{};

		/**
		 * @brief Given two pdb files containing each a model of the same protein,
		 * call TMScore binary to superimpose the first over the second. The superimposed (rotated/traslated)
		 * model is then loaded in a ProteinModel using the double pointer provided.
		 * @param modelFile (string) full path to model#1 file
		 * @param nativeFile (string) full path to model#2 file
		 * @param imposedModel(ProteinModel**) double pointer of type ProteinModel, as output
		 * @return (double) TMScore
		 */
		virtual double TMScore(string modelFile, string nativeFile, ProteinModel** imposedModel);

		/**
		 * @brief Given a ProteinModel call TMScore binary to superimpose two models contained in it.
		 * The superimposed (rotated/traslated) model is then loaded in a ProteinModel using the double pointer provided.
		 * @param prot(ProteinModel&) reference to ProteinModel object
		 * @param model (unsigned int) model#1 name in ProteinModel object
		 * @param native (unsigned int) model#2 name in ProteinModel object
		 * @param imposedModel(ProteinModel**) double pointer of type ProteinModel, as output
		 * @return (double) TMScore
		 */
		virtual double TMScore(ProteinModel& prot, unsigned int model, unsigned int native, ProteinModel** imposedModel);


		/**
		 * @brief Superimposition of two models through call to TMScore binary
		 * The superimposed (rotated/traslated) model is then loaded in a ProteinModel using the double pointer provided.
		 * @param prot1(ProteinModel&) reference to the fist ProteinModel
		 * @param model1 (unsigned int) model#1 name in the first ProteinModel object
		 * @param prot2(ProteinModel&) reference to the second ProteinModel
		 * @param model2 (unsigned int) model#2 name in the second ProteinModel object
		 * @param imposedModel(ProteinModel**) double pointer of type ProteinModel, as output
		 * @return (double) TMScore
		 */
		virtual double TMScore(ProteinModel& prot1, unsigned int model1, ProteinModel& prot2, unsigned int model2, ProteinModel** imposedModel);


		/**
		 * Set verbosity.
		 * Verbosity is applied on PdbLoader instances.
		 * @param v (bool) new verbosity option
		 */
		void setVerbose(bool v){
			verbose = v;
		}

		/**
		 * Default deconstructor
		 */
		virtual ~TMScoreBin(){};


	private:
		std::string binary;
		std::string tmp;
		bool verbose;
	};

}}



#endif /* MOBI_SOURCES_TMSCOREBIN_H_ */
