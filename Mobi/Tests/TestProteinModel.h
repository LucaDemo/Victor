/*
 * TestProteinModel.h
 *
 *  Created on: 19/giu/2015
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
#include <ProteinModel.h>

using namespace std;
using namespace Victor::Mobi;

class TestProteinModel : public CppUnit::TestFixture {

public:

	TestProteinModel(){}

	virtual ~TestProteinModel(){}

	static CppUnit::Test *suite() {
	        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestProteinModel");

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestProteinModel>("Test2 - Load and Get models",
	                &TestProteinModel::testLoadAndGet));

	        return suiteOfTests;
	    }

	void setUp(){}

	void tearDown(){}

protected:

	void testLoadAndGet(){
		cout << endl << ">>>\tTestProteinModel >>> test Load and Get models" << endl;
		string path = getenv("VICTOR_ROOT");
		string inputFile = path + "Mobi/Tests/data/1AB2_input.pdb";
		ifstream inFile(inputFile.c_str());
		if (!inFile)
			ERROR("Input file not found.", exception);

		PdbLoader pl(inFile);
		//pl.setNoVerbose();
	    ProteinModel prot = ProteinModel();
	    vector<unsigned int> models = vector<unsigned int>();
	    models.push_back(1);
	    models.push_back(4);
	    models.push_back(8);
	    cout << "Loading models 1, 4 and 8" << endl;
	    prot.load(pl,0,models);
	    cout << "Getting models back" << endl;
	    CPPUNIT_ASSERT(prot.size() == 3);
	    Spacer* m1 = prot.getModel(0);
	    Spacer* m2 = prot.getModel(1);
	    Spacer* m3 = prot.getModel(2);
	    CPPUNIT_ASSERT(m1->getAmino(0).getType() == "GLY");
	    CPPUNIT_ASSERT(m2->getAmino(1).getType() == "SER");
	    CPPUNIT_ASSERT(m3->getAmino(2).getType() == "GLY");
	    cout << "All fine!" << endl;
	}
};
