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

using namespace Victor;
using namespace Victor::Biopool;
using namespace Victor::Mobi;

void sShowHelp() {
    cout << "Help message goes here...\n\n";
}

int main(int argc, char* argv[]) {

	//If -h specified, show help message
    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }

    string inputFile, outputFile, chainID;
    unsigned int modelNum;
    bool chi, all;

    getArg("i", inputFile, argc, argv, "!");
    getArg("o", outputFile, argc, argv, "!");
    getArg("c", chainID, argc, argv, "!");
    getArg("m", modelNum, argc, argv, 999);
    all = getArg("-all", argc, argv);
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
    // Check chain args
    if ((chainID != "!") && all) {
        ERROR("You can use --all or -c, not both", error);
    }
    // User selected chain
    if (chainID != "!") {
        if (chainID.size() > 1)
            ERROR("You can choose only 1 chain", error);
        pl.setChain(chainID[0]);
    }// All chains
    else if (all) {
        pl.setAllChains();
    }// First chain
    else {
        pl.setChain(pl.getAllChains()[0]);
    }


    // Load the protein object

    // Open the proper output stream (file or stdout)
    std::ofstream fout;


    ProteinModel prot1;
    PdbSaver ps(cout);
    //ps.saveProtein(prot1);
    prot1.load(pl);
    cout << "sizeProtein = " << prot1.sizeProtein() << endl;

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
    tm.tms(prot1, 1, 2, sss);





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
