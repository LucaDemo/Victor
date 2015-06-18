/*
 * TestSDResults.h
 *
 *  Created on: 17/giu/2015
 *      Author: luca
 */

#ifndef MOBI_TESTS_TESTRESULTS_H_
#define MOBI_TESTS_TESTRESULTS_H_


#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <Spacer.h>
#include <PdbLoader.h>
#include <Results.h>

using namespace std;
using namespace Victor::Mobi;

class TestResults : public CppUnit::TestFixture {

public:

	TestResults(){}

	virtual ~TestResults(){}

	static CppUnit::Test *suite() {
	        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestResults");

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestResults>("Test1 - Populate the results map",
	                &TestResults::testSDPopulation));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestResults>("Test2 - Calculate scaled distance mean",
	                &TestResults::testSDMean));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestResults>("Test3 - Calculate scaled distance standard deviation",
	                &TestResults::testSDSD));

	        return suiteOfTests;
	    }

	void setUp(){}

	void tearDown(){}

protected:

	void testSDPopulation(){
		Results sdr = Results();
		fillWith3Results10(sdr);
		CPPUNIT_ASSERT(sdr.size() == 3);
		CPPUNIT_ASSERT(sdr.modelsSize() == 10);
	}

	void testSDMean(){
		Results sdr = Results();
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

	void fillWith3Results10(Results& sdr){
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



#endif /* MOBI_TESTS_TESTRESULTS_H_ */
