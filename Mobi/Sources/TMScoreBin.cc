/*
 * TMScoreBin.cc
 *
 *  Created on: 01/giu/2015
 *      Author: luca
 */

#include <TMScoreBin.h>
#include <PdbSaver.h>

using namespace Victor::Mobi;
using namespace Victor::Biopool;

double TMScoreBin::tms(Spacer& model, Spacer& native, Spacer* imposedModel){

	//Save spacers in pdb files
	std::ofstream fout;
	PdbSaver ps(fout);
	fout.open("pdb_mod.tmp");
	ps.saveSpacer(model);
	ps.endFile();
	fout.close();
	fout.open("pdb_nat.tmp");
	ps.saveSpacer(native);
	ps.endFile();
	fout.close();

	//Call TMScore binary
	return 0;

}


