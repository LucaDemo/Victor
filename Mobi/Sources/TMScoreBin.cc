/*
 * TMScoreBin.cc
 *
 *  Created on: 01/giu/2015
 *      Author: luca
 */

#include <TMScoreBin.h>
#include <PdbSaver.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <Debug.h>

using namespace Victor::Mobi;
using namespace Victor::Biopool;
using namespace std;

string TMTMP_IN1 = "tmin1.pdb.tmp";
string TMTMP_IN2 = "tmin2.pdb.tmp";
string TMTMP_OUT = "tmout.pdb.tmp";

double TMScoreBin::tms(ProteinModel& prot, unsigned int model, unsigned int native, Spacer& imposedModel){

	std::stringstream sstm;
	sstm << "TMScore between models " << model << " and " << native;
	DEBUG_MSG(sstm.str());
	//Save spacers in pdb files
	std::ofstream fout;
	PdbSaver ps(fout);
	fout.open(TMTMP_IN1.c_str());
	ps.saveSpacer(prot.getModel(model));
	ps.endFile();
	fout.close();
	fout.open(TMTMP_IN2.c_str());
	ps.saveSpacer(prot.getModel(native));
	ps.endFile();
	fout.close();

	//Call TMScore binary
	return tms("./" + TMTMP_IN1,"./" + TMTMP_IN2, imposedModel);
}

double TMScoreBin::tms(string modelFile, string nativeFile, Spacer& imposedModel){
	if (access(modelFile.c_str(), R_OK) == 0 && access(nativeFile.c_str(), R_OK) == 0){
		if(access(binary.c_str(), X_OK) == 0){
			pid_t pid;
			int status;
			pid = fork();
			if (pid < 0){
				ERROR("Unable to fork child process",error);
			} else {
				if (pid == 0){
					string command = binary + " " + modelFile + " " + nativeFile + " -o " + TMTMP_OUT;
					DEBUG_MSG("Exec-ing: " + command);
					if (execl("./TMScore",modelFile.c_str(),nativeFile.c_str(), NULL))
						ERROR("Unable to exec",error);
				} else{
					while (wait(&status) != pid);
					DEBUG_MSG("Complete.");
				}
			}
		}
		else
			ERROR("No access to " + binary + " binary!",exception);
	}
	else
		ERROR("No access to pdb files " + modelFile + " or " + nativeFile, exception);

	//Call TMScore binary
	return 0;
	return 0;
}


