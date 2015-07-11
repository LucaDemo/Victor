
#ifndef MOBI_SOURCES_VECTORCOLLECTION_H_
#define MOBI_SOURCES_VECTORCOLLECTION_H_

#include <Spacer.h>
#include <AtomCode.h>
#include <Atom.h>
#include <AtomCode.h>
#include <utility>
#include <map>
#include <cmath>
using namespace Victor::Biopool;

//extern double psi;

namespace Victor { namespace Mobi {

/**
 * @brief Class used to manage collections of Mobility results.
 * We use this class to process vectors of doubles resulting from mobility operations,
 * like scaled distance values. This class provides also methods to calculate means and std deviations
 * for corresponding (same position) values in different vectors.\n
 * \b Attention : all vectors must have the same dimension, the first value added to the collection
 * determines the length of all vectors.
 */
class VectorCollection{

public:
	/**
	 * Default constructor
	 */
	VectorCollection(){
		this->results = new std::map<int,std::vector<double> >();
	}

	/**
	 * Destructor
	 */
	~VectorCollection(){
		results->clear();
		delete results;
	}

	/**
	 * Clear the collection deleting all elements in it
	 */
	void clear(){
		results->clear();
	}

	/**
	 * Add values (for example scaled distance vector) to this collection
	 * @param id (int) values id
	 * @param result (vector<double>&) vector of values
	 */
	void addValue(int id, std::vector<double>& result){
		if (size() > 0)
			if (result.size() != vectorsSize())
				ERROR("Trying to add a result of non compatible size",exception);
		results->insert(std::make_pair(id,result));
		result.size();
	}

	/**
	 * Add values (for example scaled distance vector) to this collection
	 * @param m1 (int) first model id
	 * @param m2 (int) second model id
	 * @param result (vector<double>&) vector of values
	 */
	void addValue(int m1, int m2, std::vector<double>& result){
		if (size() > 0)
			if (result.size() != vectorsSize())
				ERROR("Trying to add a result of non compatible size",exception);
		results->insert(std::make_pair(VectorCollection::id(m1,m2),result));
	}

	/**
	 * Return values given the id
	 * @param id (int) id to search
	 * @return (vector<double>&) values, if found
	 */
	vector<double> getValue(int id){
		std::map<int,std::vector<double> >::const_iterator it = this->results->find(id);
		if (it != this->results->end())
			return it->second;
		else
			ERROR("Unable to find values with given id",exception);
			return vector<double>();	//never reached, just to hide warnings
	}

	/**
	 * Returns a list of all the model partecipating in this collection
	 * @return (vector<int>) models
	 */
	vector<int> getModels(){
		vector<int> models;
		for (std::map<int,std::vector<double> >::const_iterator it = this->results->begin(); it != this->results->end(); ++it){
			if (std::find(models.begin(), models.end(), VectorCollection::modelsFromID(it->first)[0]) == models.end())
				models.push_back(VectorCollection::modelsFromID(it->first)[0]);
			if (std::find(models.begin(), models.end(), VectorCollection::modelsFromID(it->first)[1]) == models.end())
							models.push_back(VectorCollection::modelsFromID(it->first)[1]);
		}
		return models;
	}

	/**
	 * Returns all records relative to the specified model
	 * @param model (int) the model to search
	 * @return (VectorCollection) subcollection relative to the model
	 */
	VectorCollection getValuesByModel(int model){
		VectorCollection outVC;
		std::map<int,std::vector<double> >::iterator it = this->results->begin();
		while (it != this->results->end()){
			if (it->first / delimiter == model)
				outVC.addValue(it->first, it->second);
			else
				if (it->first % delimiter == model)
					outVC.addValue(it->first, it->second);
			it++;
		}
		return outVC;
	}

	/**
	 * Returns the length of vectors in this collection
	 * @return (unsigned int) length
	 */
	unsigned int vectorsSize(){
		if (size() < 1)
			ERROR("Cannot get model lenght, since there are no models in this SDResult object",exception);
		return this->results->begin()->second.size();
	}

	/**
	 * Returns the size of this collection
	 * @return (unsigned int) size
	 */
	unsigned int size(){
		return this->results->size();
	}

	/**
	 * Read-Only Iterator to this collection
	 * @return (map::<int,std::vector<double>>::iterator) iterator
	 */
	std::map<int,std::vector<double> >::const_iterator iterator(){
		return this->results->begin();
	}

	/**
	 * Calculates the means of values in same position in the vectors
	 * @return (vector<double>) means
	 */
	vector<double> mean(){
		vector<double> mean = vector<double>(this->vectorsSize(),0.0);
		std::map<int,std::vector<double> >::const_iterator it;
		for (it = this->results->begin(); it != this->results->end(); ++it)
			for (unsigned int a = 0; a < this->vectorsSize(); a++)
				mean[a] += it->second[a];
		for (unsigned int a = 0; a < this->vectorsSize(); a++)
			mean[a] /= this->size();
		return mean;
	}

	/**
	 * Calculates the standard deviation of values in same position in the vectors
	 * @return (vector<double>) standard deviations
	 */
	vector<double> stdDev(){
		vector<double> mean = this->mean();
		vector<double> sd = vector<double>(this->vectorsSize(),0.0);
		std::map<int,std::vector<double> >::const_iterator it;
		for (it = this->results->begin(); it != this->results->end(); ++it)	//foreach distance record
			for (unsigned int a = 0; a < this->vectorsSize(); a++)	//foreach residue
				sd[a] += pow(it->second[a] - mean[a],2);	//cumulate the 2pow of distance minus mean
		for (unsigned int a = 0; a < this->vectorsSize(); a++)
			sd[a] = sqrt(sd[a] / this->size());
		return sd;
	}

	/**
	 * Calculate rmsd for all distance records inside this VectorCollection and return their average
	 * @return (double) mean rmsd
	 */
	double rmsd(){
		double rmsd = 0;
		double singleRmsd;
		std::map<int,std::vector<double> >::const_iterator it;
		for (it = this->results->begin(); it != this->results->end(); ++it){	//foreach distance record
			singleRmsd = 0;
			for (unsigned int i = 0; i < this->vectorsSize(); i++)	//foreach residue
				singleRmsd += pow(it->second[i],2);	//cumulate the 2pow of distance
			singleRmsd = sqrt(singleRmsd / this->vectorsSize());
			rmsd += singleRmsd;
		}
		return (rmsd / results->size());
	}

	/**
	 * Short output... debug purpouses
	 */
	void print(){
		std::map<int,std::vector<double> >::iterator it;
		for (it = this->results->begin(); it != this->results->end(); ++it)
			cout << it->first << " => [" << it->second[0] << "," << it->second[1] << ", ...]" << endl;
	}

	static const int delimiter = 1000;
	static int id(int m1, int m2){
		return m1*delimiter + m2;
	}
	static vector<int> modelsFromID(int id){
		vector<int> models(2);
		models[0] = id / delimiter;
		models[1] = id % delimiter;
		return models;
	}




protected:
	std::map<int,std::vector<double> > *results;
};


}}	//namespaces

#endif /* MOBI_SOURCES_VECTORCOLLECTION_H_ */
