/*
 * TestMobi.cc
 *
 *  Created on: 17/giu/2015
 *      Author: luca
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include "TestVectorCollection.h"
using namespace std;


int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites:" << endl;
	runner.addTest(TestVectorCollection::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();

	return 0;
}
