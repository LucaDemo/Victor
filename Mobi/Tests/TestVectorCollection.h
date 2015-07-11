/*
 * TestSDResults.h
 *
 *  Created on: 17/giu/2015
 *      Author: luca
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <Spacer.h>
#include <PdbLoader.h>
#include "../Sources/VectorCollection.h"

using namespace std;
using namespace Victor::Mobi;

class TestVectorCollection : public CppUnit::TestFixture {

public:

	TestVectorCollection(){}

	virtual ~TestVectorCollection(){}

	static CppUnit::Test *suite() {
	        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestResults");

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestVectorCollection>("Test1 - Populate the map",
	                &TestVectorCollection::testPopulation));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestVectorCollection>("Test2 - Calculate means",
	                &TestVectorCollection::testMeanSD));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestVectorCollection>("Test3 - Calculate standard deviation",
	                &TestVectorCollection::testSD));

	        return suiteOfTests;
	    }

	void setUp(){}

	void tearDown(){}

protected:

	void testPopulation(){
		VectorCollection sdr = VectorCollection();
		cout << endl << ">>>\tTestVectorCollection >>> test Population..." << endl;
		fillWith5Results10(sdr);
		CPPUNIT_ASSERT(sdr.size() == 5);
		CPPUNIT_ASSERT(sdr.vectorsSize() == 10);
		std::map<int,std::vector<double> >::const_iterator it = sdr.iterator();
		CPPUNIT_ASSERT(it->second[5] == 5);
		it++;
		CPPUNIT_ASSERT(it->second[5] == 6);
		it++;
		CPPUNIT_ASSERT(it->second[5] == 7);

		vector<int> ids = sdr.getModels();
		CPPUNIT_ASSERT(ids.size() == 6);
	}

	void testMeanSD(){
		cout << endl << ">>>\tTestVectorCollection >>> test Mean..." << endl;
		VectorCollection sdr = VectorCollection();
		fillWith5Results10(sdr);
		vector<double> mean = sdr.mean();
		for (unsigned int i = 0; i < 10; i++)
			CPPUNIT_ASSERT(mean[i] == i+2);
	}

	void testSD(){
		cout << endl << ">>>\tTestVectorCollection >>> test deviation calculation..." << endl;
		VectorCollection sdr = VectorCollection();
		fillWith5Results10(sdr);
		vector<double> sd = sdr.stdDev();
		for (unsigned int i = 0; i < 10; i++)
			CPPUNIT_ASSERT(sd[i] = sqrt(2));
	}

private:

	void fillWith5Results10(VectorCollection& sdr){
		unsigned int modelLen = 10;
		vector<double> *v1 = new vector<double>(modelLen);
		vector<double> *v2 = new vector<double>(modelLen);
		vector<double> *v3 = new vector<double>(modelLen);
		vector<double> *v4 = new vector<double>(modelLen);
		vector<double> *v5 = new vector<double>(modelLen);
		for (unsigned int i = 0; i < modelLen; i++){
			(*v1)[i] = i;
			(*v2)[i] = i+1;
			(*v3)[i] = i+2;
			(*v4)[i] = i+3;
			(*v5)[i] = i+4;
		}
		sdr.addValue(1,2,*v1);
		sdr.addValue(2,3,*v2);
		sdr.addValue(3,4,*v3);
		sdr.addValue(4,4,*v4);
		sdr.addValue(5,7,*v5);
	}
};

