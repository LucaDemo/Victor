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
//#include <Debug.h>

using namespace Victor::Mobi;
using namespace Victor::Biopool;
//using namespace std;

//char* TMTMP_IN1 = "tmin1.pdb.tmp";
//char* TMTMP_IN2 = "tmin2.pdb.tmp";
//char* TMTMP_OUT = "tmout.pdb.tmp";

double TMScoreBin::tms(ProteinModel& prot, unsigned int model, unsigned int native, Spacer& imposedModel){

	//Save spacers in pdb files
	std::ofstream fout;
	PdbSaver ps(fout);
	fout.open("");
	ps.saveSpacer(prot.getModel(model));
	ps.endFile();
	fout.close();
	fout.open("");
	ps.saveSpacer(prot.getModel(native));
	ps.endFile();
	fout.close();

	//Call TMScore binary
	return tms("", "", imposedModel);
}

double TMScoreBin::tms(const char* modelFile, const char* nativeFile, Spacer& imposedModel){
	if (access(modelFile, R_OK) && access(nativeFile, R_OK)){
		if(access(binary.c_str(), X_OK)){
			pid_t pid;
			int status;
			pid = fork();
			if (pid < 0){
				//ERROR("Unable to fork child process",error);
			} else {
				if (pid == 0){
					//DEBUG_MSG("Exec-ing: " + binary.c_str() + " " + modelFile + " " + nativeFile + " -o " + *TMTMP_OUT);
					if (execlp(binary.c_str(), modelFile, nativeFile, "-o ") < 0)
						cout << "error";
						//ERROR("Unable to exec",error);
				} else{
					while (wait(&status) != pid);
					//DEBUG_MSG("Complete.");
				}
			}
		}
	}

	//Call TMScore binary
	return 0;
	return 0;
}


