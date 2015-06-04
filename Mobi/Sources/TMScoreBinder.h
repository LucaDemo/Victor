/*
 * TMScoreBinder.h
 *
 *  Created on: 01/giu/2015
 *      Author: luca
 */

#ifndef MOBI_SOURCES_TMSCOREBINDER_H_
#define MOBI_SOURCES_TMSCOREBINDER_H_

#include<Spacer.h>
#include<Protein.h>

using namespace Victor::Biopool;

namespace Victor { namespace Mobi {
	/**
	 * @brief base class that represent a binder to TMScore functionalities.
	 */
	class TMScoreBinder{
	public:
		/**
		 * @brief Calculates TMScore with given models
		 * @param prot (ProteinModel&) protein object
		 * @param model (unsigned int) model name (the one to be super-imposed)
		 * @param native (unsigned int) native model name (the one on which we super-impose)
		 * @newMode (Spacer*) target model aligned to reference
		 * @return TMScore of the two models
		 */
		virtual double tms (ProteinModel& prot, unsigned int model, unsigned int native, Spacer& imposedModel){return 0;};
		virtual ~TMScoreBinder(){};
	};

}}



#endif /* MOBI_SOURCES_TMSCOREBINDER_H_ */
