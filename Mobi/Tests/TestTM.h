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
#include <TMScoreBin.h>

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

		PdbLoader pl(inFile);
		pl.setNoVerbose();
	    ProteinModel prot;
	    vector<unsigned int> models = vector<unsigned int>();
	    models.push_back(1);
	    models.push_back(2);
	    cout << "Loading models 1 and 2" << endl;
	    prot.load(pl,0,models);
	    cout << "TM Score binary call" << endl;
	    TMScoreBin tmsb = TMScoreBin(TMDir+"TMScore", TMDir);
	    Spacer* sp;
//	    double score = tmsb.tms(prot,0,1,&sp);
//	    cout << "Score is " << score << endl;
//	    CPPUNIT_ASSERT(score > 0 && score < 1);

	    ifstream outFile((TMDir + TMTMP_OUT + "_atm").c_str());
	    //while (outFile)
		//	cout << readLine(outFile);


	    PdbLoader pl2(outFile);
	    cout << pl2.getMaxModels() << " models" << endl;
	    pl2.setModel(2);
	    pl2.checkModel();
	    ProteinModel prot2;
	    models.clear();
	    models.push_back(2);
	    prot2.load(pl2);

//	    CPPUNIT_ASSERT(prot.size() == 3);
//	    Spacer m1 = prot.getModel(0);
//	    Spacer m2 = prot.getModel(1);
//	    Spacer m3 = prot.getModel(2);
//	    CPPUNIT_ASSERT(m1.getAmino(0).getType() == "GLY");
//	    CPPUNIT_ASSERT(m2.getAmino(1).getType() == "SER");
//	    CPPUNIT_ASSERT(m3.getAmino(2).getType() == "GLY");
//	    cout << "All fine!" << endl;
	}
};


