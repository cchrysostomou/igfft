#define _CRT_SECURE_NO_DEPRECATE

#ifndef __ReadFiles_h__
#define __ReadFiles_h__

#include <iostream>
#include <algorithm>
#include <fstream> //for ifstream
#include <cstring> //for strlen
#include <cstdio>
#include <vector>
#include <typeinfo> //for typeid
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string>
#include <sstream> //to use istringstream and getline
#include <functional>
#include <cctype>
#include <locale>

#include "StructVars.h"

namespace readfiles{
	//just extract one sequence at a time from a file
	structvars::FASTQinfo readIndividualSequence(std::ifstream &);
	structvars::FASTAinfo readIndividualSequenceToFASTA(std::ifstream &);
	structvars::FASTQinfo readIndividualFASTQSequence(std::ifstream &);
	structvars::FASTQinfo readIndividualTabSequence(std::ifstream &);

	//store all of the sequences from a file into a vector
	structvars::Fileinfo readFASTAToVector(const char*, std::vector<structvars::FASTQinfo>&);

	//store all of the sequences from a file into a vector
	structvars::Fileinfo readFASTQToVector(const char*, std::vector<structvars::FASTQinfo>&,char filetype='Q');

	//store all sequences from tab delimted file as vector
	structvars::Fileinfo readTABToVector(const char*, std::vector<structvars::FASTQinfo>&);

	//scan through a fasta file, fastq file, or a tab delimited file and return the number of sequences and the lenght of the maximum sequence and length of the minimum esquence.
	structvars::Fileinfo readSeqLengths(const char*, std::string,bool ignoreHeader=true);

	//convert a fasta file to a tab delimited file
	structvars::Fileinfo convertFASTAtoTAB(std::string &);

	//convert a fastq file to a tab delimited file
	structvars::Fileinfo convertFASTQtoTAB(std::string &);

	//take a fasta file and make it into a tab delimited file
	void WriteFASTASeqToTab(std::ofstream &, structvars::FASTAinfo);
	//take a fasta file and make it into a tab delimited file
	void WriteFASTASeqToTab(std::ofstream &, structvars::FASTQinfo);

	//take a fastq file and make it into a tab delimited file
	void WriteFASTQSeqToTab(std::ofstream &, structvars::FASTQinfo);

	//to parse the parwisealignment struct string....see cpp file:
	structvars::PairwiseClusterAlignment ParsePairwiseAlignmentString(std::string);


	//search for a substring
	void FindAbRegion(structvars::Abregion &, std::string, int &);

	//remove trailing spaces
	void RemoveTrailingSpaces(std::string & s);

	//update the vector containing binned sequence lengths
	//we will use this vector to determine the length distributions (i.e. where the majority of the sequences are sized)
	inline void AddSeqLenToVecDistribution(int len, std::vector<std::vector<int>> & lengthDist){
		bool found = false;		
		int base = int(log2(len)),minb=pow(2,base);
		double precision = 0.25,inverse=1/precision, binNum = precision*1.0*minb*floor(inverse*(((len*1.0) / (minb)) - 1));
		int index = minb + binNum, upper = index+minb/inverse - 1;
		for (int i = 0; i < lengthDist.size(); i++){
			if (lengthDist[i][0] == index){
				lengthDist[i][3]++;
				lengthDist[i][2] += len;
				found = true;
				break;
			}
		}
		if (!found)
			lengthDist.push_back(std::vector<int>{ index, upper , len, 1 });
	}

}
#endif

