/*
 * TestSDResults.h
 *
 *  Created on: 17/giu/2015
 *      Author: luca
 */

#ifndef MOBI_TESTS_TESTSDRESULTS_H_
#define MOBI_TESTS_TESTSDRESULTS_H_


#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <Spacer.h>
#include <PdbLoader.h>
#include <SDResults.h>

using namespace std;
using namespace Victor::Mobi;

class TestSDResults : public CppUnit::TestFixture {

public:

	TestSDResults(){}

	~TestSDResults(){}

	static CppUnit::Test *suite() {
	        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestSDResults");

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestSDResults>("Test1 - Populate the results map",
	                &TestSDResults::testSDPopulation));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestSDResults>("Test2 - Calculate scaled distance mean",
	                &TestSDResults::testSDMean));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestSDResults>("Test3 - Calculate scaled distance standard deviation",
	                &TestSDResults::testSDSD));

	        return suiteOfTests;
	    }

	void setUp(){}

	void tearDown(){}

protected:

	void testSDPopulation(){
		SDResults sdr = SDResults();
		fillWith3Results10(sdr);
		CPPUNIT_ASSERT(sdr.size() == 3);
		CPPUNIT_ASSERT(sdr.modelsSize() == 10);
	}

	void testSDMean(){
		SDResults sdr = SDResults();
		fillWith3Results10(sdr);
		vector<double> *mean = new vector<double>(10);
		sdr.meanSD(mean);
		for (unsigned int i = 0; i < 10; i++)
			CPPUNIT_ASSERT((*mean)[i] == i+3);
	}

	void testSDSD(){
		CPPUNIT_ASSERT(true);
	}

private:

	void fillWith3Results10(SDResults& sdr){
		unsigned int modelLen = 10;
		vector<double> *v1 = new vector<double>(modelLen);
		vector<double> *v2 = new vector<double>(modelLen);
		vector<double> *v3 = new vector<double>(modelLen);
		for (unsigned int i = 0; i < modelLen; i++){
			(*v1)[i] = i;
			(*v2)[i] = i+1;
			(*v3)[i] = i+2;
		}
		sdr.addResult(1,*v1);
		sdr.addResult(2,*v2);
		sdr.addResult(3,*v3);
	}
};



#endif /* MOBI_TESTS_TESTSDRESULTS_H_ */
