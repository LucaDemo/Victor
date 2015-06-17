/*
 * TMScoreBinder.h
 *
 *  Created on: 01/giu/2015
 *      Author: luca
 */

#ifndef MOBI_SOURCES_TMSCOREBINDER_H_
#define MOBI_SOURCES_TMSCOREBINDER_H_

#include<Spacer.h>
#include<ProteinModel.h>

using namespace Victor::Biopool;
using namespace Victor::Mobi;

namespace Victor { namespace Mobi {
	/**
	 * @brief base class that represent a binder to TMScore functionalities.
	 */
	class TMScoreBinder{
	public:
		/**
		 * @brief Calculates TMScore with given models
		 * @param prot (ProteinModel&) protein object containing the models
		 * @param model (unsigned int) target model (the one to be super-imposed)
		 * @param native (unsigned int) native model (the one on which we super-impose)
		 * @param imposedModel (Spacer*) copy of target model after TMScore super-imposition over native model
		 * @return TMScore of the two models
		 */
		virtual double tms (ProteinModel& prot, unsigned int model, unsigned int native, Spacer& imposedModel){return 0;};

		virtual double tms(const char* modelFile, const char* nativeFile, Spacer& imposedModel);
		virtual ~TMScoreBinder(){}
	};

}}



#endif /* MOBI_SOURCES_TMSCOREBINDER_H_ */
