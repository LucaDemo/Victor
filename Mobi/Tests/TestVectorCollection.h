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
		cout << endl << ">>>\tTestVectorCollection >>> test Population\nPopulating collection..." << endl;
		fillWith3Results10(sdr);
		CPPUNIT_ASSERT(sdr.size() == 3);
		CPPUNIT_ASSERT(sdr.vectorsSize() == 10);
		cout << "Correctly populated!" << endl;
		std::map<int,std::vector<double> >::const_iterator it = sdr.iterator();
		CPPUNIT_ASSERT(it->second[5] == 5);
		it++;
		CPPUNIT_ASSERT(it->second[5] == 6);
		it++;
		CPPUNIT_ASSERT(it->second[5] == 7);
		cout << "Correctly retrieved!" << endl;
	}

	void testMeanSD(){
		cout << endl << ">>>\tTestVectorCollection >>> test Mean\nPopulating collection..." << endl;
		VectorCollection sdr = VectorCollection();
		fillWith3Results10(sdr);
		vector<double> mean = sdr.mean();
		for (unsigned int i = 0; i < 10; i++)
			CPPUNIT_ASSERT(mean[i] == i+1);
		cout << "Mean calculation correct!" << endl;
	}

	void testSD(){
		cout << endl << ">>>\tTestVectorCollection >>> test deviation calculation\nPopulating collection..." << endl;
		VectorCollection sdr = VectorCollection();
		fillWith3Results10(sdr);
		vector<double> sd = sdr.stdDev();
		for (unsigned int i = 0; i < 10; i++){
			cout << sd[i] << endl;
			CPPUNIT_ASSERT(sd[i] > 0.816496 && sd[i] < 0.816498);
		}

	}

private:

	void fillWith3Results10(VectorCollection& sdr){
		unsigned int modelLen = 10;
		vector<double> *v1 = new vector<double>(modelLen);
		vector<double> *v2 = new vector<double>(modelLen);
		vector<double> *v3 = new vector<double>(modelLen);
		for (unsigned int i = 0; i < modelLen; i++){
			(*v1)[i] = i;
			(*v2)[i] = i+1;
			(*v3)[i] = i+2;
		}
		sdr.addValue(1,*v1);
		sdr.addValue(2,*v2);
		sdr.addValue(3,*v3);
	}
};

