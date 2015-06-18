/*
 * TestMobi.cc
 *
 *  Created on: 17/giu/2015
 *      Author: luca
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include <TestSDResults.h>
using namespace std;


int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites:" << endl;
	runner.addTest(TestSDResults::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();

	return 0;
}
