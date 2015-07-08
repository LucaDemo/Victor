/* 	This file is part of Victor

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
 */

/**
 * @file TMScoreBin.cc
 * @author Luca Demo
 * @date Jun 2015
 * @version 0.1
 */
#include <TMScoreBin.h>
#include <PdbSaver.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <Debug.h>
#include <stdio.h>
#include <string>
#include <iostream>


using namespace Victor::Mobi;
using namespace Victor::Biopool;
using namespace std;

const string TMTMP_IN1 = "tmin1.pdb.tmp";
const string TMTMP_IN2 = "tmin2.pdb.tmp";
const string TMTMP_OUT = "tmout.pdb.tmp";

/**
 * @brief the TMScore output is not pdb conformant. This static method read the output and fix it in a memory buffer.
 * Then loads a ProteinModel with the superimposed model only.
 * @param pdbFile (string) full path to TMScore output
 * @param imposedModel (ProteinModel**) double pointer of type PRoteinModel, as output
 */
void spacerFromTMOutput(string pdbFile, ProteinModel** imposedModel);


void TMScoreBin::TMImpose(ProteinModel& prot1, unsigned int model1, ProteinModel& prot2, unsigned int model2, ProteinModel** imposedModel){
	if (verbose)
		cout << "TMScore between models " << model1 << " and " << model2;

	//Save spacers in pdb files
	std::ofstream fout;
	PdbSaver ps(fout);
	fout.open((tmp + TMTMP_IN1).c_str());
	ps.saveSpacer(*(prot1.getModel(model1)));
	ps.endFile();
	fout.close();
	fout.open((tmp + TMTMP_IN2).c_str());
	ps.saveSpacer(*(prot2.getModel(model2)));
	ps.endFile();
	fout.close();
	//Call TMScore binary
	return TMImpose((tmp + TMTMP_IN1), (tmp + TMTMP_IN2), imposedModel);
}

void TMScoreBin::TMImpose(ProteinModel& prot, unsigned int model1, unsigned int model2, ProteinModel** imposedModel){
	return TMImpose(prot,model1,prot,model2,imposedModel);
}

void TMScoreBin::TMImpose(string modelFile, string nativeFile, ProteinModel** imposedModel){
	double score = -1;
	if (access(modelFile.c_str(), R_OK) == 0 && access(nativeFile.c_str(), R_OK) == 0){
		if(access(binary.c_str(), X_OK) == 0){
			pid_t pid;
			//Pipeing TMS output
			int pipefd[2];
			pipe(pipefd);
			//Forking child for TMS
			pid = fork();
			if (pid < 0){
				ERROR("Unable to fork child process",error);
			} else {
				if (pid == 0){
					close(pipefd[0]);
					dup2(pipefd[1],1);	//stdout
					dup2(pipefd[1],2);	//stderr
					close(pipefd[1]);
					if (execl(binary.c_str(),binary.c_str(), modelFile.c_str(),nativeFile.c_str(), "-o" , (tmp + TMTMP_OUT).c_str(),NULL))
						ERROR("Unable to exec",error);
				} else{
					//Read output from child pipe
					char buffer[2048];
					if (read(pipefd[0],buffer, sizeof(buffer)) > 0){
						char * tok;
						tok = strtok(buffer,"\n\r");
						while (tok != NULL){
							//cout << string(tok) << endl;
							if (string(tok).find("TM-score") == 0)
								if (string(tok).find("=") > 0)
									score = std::strtod(string(tok).substr(string(tok).find("=")+2,6).c_str(),NULL);
							tok = strtok(NULL,"\n\r");
						}
					}
					if (score < 0)
						ERROR("TM-score failed or non recognised output", exception)
				}
			}
		}
		else
			ERROR("No access to " + binary + "  binary!",exception);
	}
	else
		ERROR("No access to pdb files " + modelFile + " or " + nativeFile, exception);

	return spacerFromTMOutput(tmp + TMTMP_OUT + "_atm", imposedModel);
}

void spacerFromTMOutput(string pdbFile, ProteinModel** imposedModel){
	if (access(pdbFile.c_str(), R_OK) != 0)
		ERROR("Cannot read pdb file to fix",exception);

	ifstream inFile(pdbFile.c_str());
	stringstream buffer;
	string line;
	while (inFile){
		line = readLine(inFile);
		if (line.substr(0,16) == "REMARK  TM-score"){
			buffer << line << endl;
			if (line.substr(6,8) == "TM-score")
				buffer << "MODEL        1" << endl;
		}else if (line.substr(0,6) == "ATOM  "){
			buffer << line << endl;
		}else if (line.substr(0,3) == "TER"){
			buffer << line << endl;
			buffer << "ENDMDL" << endl;
			break;
		}
	}
	buffer.clear();
	buffer.seekg(0);
	//Load from memory buffer
	PdbLoader pl(buffer);
	pl.setNoVerbose();
	//Load and return protein object
	*imposedModel = new ProteinModel();
	pl.setModel(1);
	cout.flush();
	pl.loadProtein(**imposedModel);
}

