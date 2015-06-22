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
	class TMScoreBin{
	public:

		/**
		 * @brief Setup.
		 * @param _binary (string) full path to binary TMScore file, must have execution permission
		 * @param _tmp (string) full path to temp dir, must have write permission
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
		 */
		virtual void TMImpose(string modelFile, string nativeFile, ProteinModel** imposedModel);

		/**
		 * @brief Given a ProteinModel call TMScore binary to superimpose two models contained in it.
		 * The superimposed (rotated/traslated) model is then loaded in a ProteinModel using the double pointer provided.
		 * @param prot(ProteinModel&) reference to ProteinModel object
		 * @param model (unsigned int) model#1 name in ProteinModel object
		 * @param native (unsigned int) model#2 name in ProteinModel object
		 * @param imposedModel(ProteinModel**) double pointer of type ProteinModel, as output
		 */
		virtual void TMImpose(ProteinModel& prot, unsigned int model, unsigned int native, ProteinModel** imposedModel);

		/**
		 * @brief the TMScore output is not pdb conformant. This static method read the output and fix it in a memory buffer.
		 * Then loads a ProteinModel with the superimposed model only.
		 * @param pdbFile (string) full path to TMScore output
		 * @param imposedModel (ProteinModel**) double pointer of type PRoteinModel, as output
		 */
		virtual void spacerFromTMOutput(string pdbFile, ProteinModel** imposedModel);

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
		~TMScoreBin(){};


	private:
		std::string binary;
		std::string tmp;
		bool verbose;
	};

}}



#endif /* MOBI_SOURCES_TMSCOREBIN_H_ */
