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

//extern const char* TMTMP_IN1;
//extern const char* TMTMP_IN2;
//extern const char* TMTMP_OUT;


namespace Victor { namespace Mobi {
	/**
	 * @brief TMScore functionalities through external binary.
	 */
	class TMScoreBin : public TMScoreBinder{
	public:

		/**
		 * @brief Setup.
		 * @param _binary (string) full path to binary TMScore file, must have execution permission
		 * @param _tmp (string) full path to temp dir, must have write permission
		 */
		TMScoreBin(std::string _binary = "TMScore", std::string _tmp = ".") : binary(_binary), tmp(_tmp){};


		double tms(const char* modelFile, const char* nativeFile, Spacer& imposedModel);
		double tms(ProteinModel& prot, unsigned int model, unsigned int native, Spacer& imposedModel);
		~TMScoreBin(){
			this->~TMScoreBinder();
		};

	private:
		std::string binary;
		std::string tmp;
	};

}}



#endif /* MOBI_SOURCES_TMSCOREBIN_H_ */
