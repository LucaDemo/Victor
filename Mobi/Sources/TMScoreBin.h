/*
 * TMScoreBin.h
 *
 *  Created on: 01/giu/2015
 *      Author: luca
 */

#ifndef MOBI_SOURCES_TMSCOREBIN_H_
#define MOBI_SOURCES_TMSCOREBIN_H_

#include <TMScoreBinder.h>


using namespace Victor::Biopool;
using namespace Victor::Mobi;

extern string TMTMP_IN1;
extern string TMTMP_IN2;
extern string TMTMP_OUT;


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
		TMScoreBin(std::string _binary = "TMScore", std::string _tmp = ".") : binary(_binary), tmp(_tmp){};


		virtual double tms(string modelFile, string& nativeFile, Spacer& imposedModel);
		virtual double tms(ProteinModel& prot, unsigned int model, unsigned int native, Spacer& imposedModel);
		virtual ~TMScoreBin(){};


	private:
		std::string binary;
		std::string tmp;
	};

}}



#endif /* MOBI_SOURCES_TMSCOREBIN_H_ */
