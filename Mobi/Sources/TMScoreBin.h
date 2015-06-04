/*
 * TMScoreBin.h
 *
 *  Created on: 01/giu/2015
 *      Author: luca
 */

#ifndef MOBI_SOURCES_TMSCOREBIN_H_
#define MOBI_SOURCES_TMSCOREBIN_H_

#include <TMScoreBinder.h>
#include <Spacer.h>

using namespace Victor::Biopool;
using namespace Victor::Mobi;


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
		TMScoreBin(std::string _binary, std::string _tmp) : binary(_binary), tmp(_tmp){};

		double tms(Spacer& model, Spacer& native, Spacer* imposedModel);
		~TMScoreBin(){

		};

	private:
		std::string binary;
		std::string tmp;
	};

}}



#endif /* MOBI_SOURCES_TMSCOREBIN_H_ */
