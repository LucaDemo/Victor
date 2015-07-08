/*
 * TestMobi.cc
 *
 *  Created on: 17/giu/2015
 *      Author: luca
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include <TestProteinModel.h>
#include <TestVectorCollection.h>
#include <TestTM.h>
#include <TestMobiMethods.h>
#include <MadTests.h>
using namespace std;


int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites:" << endl;
//	runner.addTest(TestVectorCollection::suite());
//	runner.addTest(TestProteinModel::suite());
//	runner.addTest(TestTM::suite());
//	runner.addTest(TestMobiMethods::suite());
	runner.addTest(MadTests::suite());
	cout<< "Running the unit tests. " << endl;

	runner.run();
	return 0;
}
