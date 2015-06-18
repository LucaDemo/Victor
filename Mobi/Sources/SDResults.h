
#ifndef MOBI_SOURCES_SDRESULTS_H_
#define MOBI_SOURCES_SDRESULTS_H_

#include <Spacer.h>
#include <AtomCode.h>
#include <Atom.h>
#include <AtomCode.h>
#include <utility>
#include <map>
using namespace Victor::Biopool;

extern const double DEF_D0;
extern const AtomCode DEF_ATOM;
//extern double psi;

namespace Victor { namespace Mobi {

void scaledDistance(vector<double>& sd, Spacer& mod, Spacer& ref, AtomCode atom = CA, double d0 = 4);
//void avgScaledDistance();

/**
 * @brief Class used to manage collections of Mobility results
 */
class SDResults{

public:
	SDResults(){
		this->results = new std::map<int,std::vector<double> >();
	}

	void addResult(int id, std::vector<double>& result){
		if (size() > 0)
			if (result.size() != modelsSize())
				ERROR("Trying to add a result of non compatible size",exception);
		results->insert(std::make_pair(id,result));
		result.size();
	}

	void getResult(vector<double>& res, int id){
		std::map<int,std::vector<double> >::const_iterator it = this->results->find(id);
		if (it != this->results->end())
			res = it->second;
	}

	unsigned int modelsSize(){
		if (size() < 1)
			ERROR("Cannot get model lenght, since there are no models in this SDResult object",exception);
		return this->results->begin()->second.size();
	}

	unsigned int size(){
		return this->results->size();
	}

	std::map<int,std::vector<double> >::const_iterator iterator(){
		return this->results->begin();
	}

	void meanSD(vector<double> *res);

private:
	std::map<int,std::vector<double> > *results;
};


}}	//namespaces

#endif /* MOBI_SOURCES_SDRESULTS_H_ */
