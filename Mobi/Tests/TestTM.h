/*
 * TestTM.h
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
#include <PdbSaver.h>
#include <TMScoreBin.h>
#include <AminoAcid.h>

using namespace std;
using namespace Victor::Mobi;

class TestTM : public CppUnit::TestFixture {

public:

	TestTM(){}

	virtual ~TestTM(){}

	static CppUnit::Test *suite() {
	        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestTM");

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestTM>("Test1 - TM Score",
	                &TestTM::testLoadAndGet));

	        return suiteOfTests;
	    }

	void setUp(){}

	void tearDown(){}

protected:

	void testLoadAndGet(){
		cout << endl << ">>>\tTestTM >>> test TM binary call" << endl;
		string path = getenv("VICTOR_ROOT");
		string inputFile = path + "Mobi/Tests/data/1AB2_input.pdb";
		string TMDir = path + "Mobi/Tests/data/";
		ifstream inFile(inputFile.c_str());
		if (!inFile)
			ERROR("Input file not found.", exception);
		//models to load
		vector<unsigned int> models = vector<unsigned int>();
	    models.push_back(3);
	    models.push_back(4);

	    //loading models
		PdbLoader pl(inFile);
		pl.setNoVerbose();
		ProteinModel prot;
	    prot.load(pl,models);

	    //calling TMScore bin
	    TMScoreBin tmsb(TMDir+"TMScore", TMDir, false);
	    ProteinModel* superImposed;
	    tmsb.TMScore(prot,0,1,&superImposed);


	    //save superimposed model to file
	    std::ofstream fout;
		PdbSaver ps(fout);
		fout.open((TMDir + "superimposed_TM-Test.pdb").c_str());
		ps.saveSpacer(*(superImposed->getModel(0)));
		ps.endFile();
		fout.close();

		CPPUNIT_ASSERT(superImposed->getModel(0)->sizeAmino() == 109);
		CPPUNIT_ASSERT(superImposed->getModel(0)->getAmino(0).getType() == "GLY");
	}
};


