/*  This file is part of Victor.

    Victor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Victor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Victor.  If not, see <http://www.gnu.org/licenses/>.

    This file is part of Mobi extension created by
 	Luca Demo  -  University of Padova, department of Science
 */
/**
 */

#include <Protein.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IoTools.h>
#include <GetArg.h>
#include <ProteinModel.h>
#include <TMScoreBin.h>
#include <Debug.h>
#include <String2Number.h>

using namespace Victor;
using namespace Victor::Biopool;
using namespace Victor::Mobi;

void sShowHelp() {
    cout << "Help message goes here...\n\n";
}

vector<std::string> splitString(std::string in, char del){
	vector<std::string> toks;
	std::stringstream ss(in);
	std::string tok;

	while(std::getline(ss, tok, del)){
		toks.push_back(tok);
	}
	return toks;
}

int main(int argc, char* argv[]) {

	DEBUG_MSG("Welcome to Mobi!");
	//If -h specified, show help message
    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }

    string inputFile, outputFile, chainID, modelList;
    vector<unsigned int> models;
    unsigned int modevelNum;

    bool chi;

    getArg("i", inputFile, argc, argv, "!");
    getArg("o", outputFile, argc, argv, "!");
    getArg("c", chainID, argc, argv, "!");
    getArg("m", modelList, argc, argv, "!");
    chi = getArg("-chi", argc, argv);

    // Check input file
    if (inputFile == "!") {
        cout << "Missing input file specification. Aborting. (-h for help)" << endl;
        return -1;
    }
    ifstream inFile(inputFile.c_str());
    if (!inFile)
        ERROR("Input file not found.", exception);

    PdbLoader pl(inFile);
    vector<Atom*> model;

    // Set PdbLoader variables
    if (!getArg("v", argc, argv)) {
        pl.setNoVerbose();
    }

    // User selected chain
    if (chainID != "!") {
        if (chainID.size() > 1)
            ERROR("You can choose only 1 chain", error);
        pl.setChain(chainID[0]);
        pl.checkAndSetChain();
    }// First chain
    else {
        pl.setChain(pl.getAllChains()[0]);
    }

    // User selected models
    if (modelList != "!"){
    	models = sToVectorOfUIntDEF(translate(modelList,',',' '));
    		for (unsigned int i = 0; i < models.size(); i++)
    			if (models[i] > pl.getMaxModels())
    				ERROR("You specified out of bound model number(s). Max for this input pdb is " + pl.getMaxModels(), exception);
    }else{
    	for (unsigned int i = 1; i <= pl.getMaxModels(); i++)
    		models.push_back(i);
    }


    // Load the protein object

    // Open the proper output stream (file or stdout)
    //std::ofstream fout;


    ProteinModel prot;
    //PdbSaver ps(cout);
    models = vector<unsigned int>();
    models.push_back(3);
    models.push_back(1);
    //models.push_back(3);
    prot.load(pl,models);


    /*
    for (unsigned int i = 1; i < 5; i++){
    	std::ostringstream ss;
    	ss << "pdbout" << i << ".pdb";
    	fout.open(ss.str().c_str());
    	ps.saveSpacer(prot1.getModel(i));
    	ps.endFile();
    	fout.close();
    }
	*/
    TMScoreBin tm("TMScore",".");
    Spacer sss;
    cout << "TM-score = " << tm.tms(prot, 0, 1, sss);


    //dtf


    /*
    ProteinModel prot2;
    fout.open("pdbout2.pdb");
	//pl.setModel(2);
	prot2.load(pl);
    //prot2.load(pl);
	cout << "sizeProtein = " << prot2.sizeProtein() << endl;
	ps.saveSpacer(prot2.getModel(0));
	ps.endFile();
	fout.close();
	*/
    return 0;
}
