#define _CRT_SECURE_NO_DEPRECATE

#ifndef __DNAFunctions_h__
#define __DNAFunctions_h__

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
#include <numeric>

#include "FFTFunctions.h"
#include "StructVars.h"

struct compute_algn_score{
	int operator() (const char &x, const char &y)const {
		if ((x != NOT_NT) && (y != NOT_NT))
			return x == y ? match : mismatch;
		else
			return 0;
	}
	typedef int first_argument_type;
	typedef int second_argument_type;
	typedef int result_type;
	int match, mismatch;
	compute_algn_score(int m, int mm){
		match = m;
		mismatch = -1 * mm;
	}
};





//%THIS FUNCTION JUST COUNTS THE NUMBER OF MATCHES betwen sequences already
/*%aligned to one another. if the number of matches is below a threshold , then we return true. if not, then false*/
//
//matchThresh = minimum number of base pair matches for boolean true
//algnPos[0] => seqStart
//algnPos[1] = > germStart
//algnPos[2] => algnLen
//algnPos[3] = > seqLen
//algnPos[4] = > germLen



namespace dnafunctions{
	

	std::string ReverseComplement(std::string &); //CACLUATE REVERSE COMPLEMENT. ALSO ANY CHARACTERS THAT ARE NOT ACTG WILL BE CONVERTED TO NULL (ASCI==0)
		
	const std::map<std::string, char> GenerateCodonMap(); //generate a codon map (codons to amion acid) for creating 
	const std::map<std::string, char> dna_to_aa_map = GenerateCodonMap();
	
	std::string TranslateSeq(const std::string & nt_sequence,int frame=1,int start =0); //convert sequences to amino acid 

	
	inline void PadSequence(std::string & dnastring, int newLength, int direction){
		int seqLength = dnastring.length();
		size_t addedLetters = newLength - seqLength;
		if (direction == FRONT)
			dnastring.insert(0, addedLetters, NOT_NT);
		else if (direction == END)
			dnastring.append(addedLetters, NOT_NT);
		return;
	}

	double AlignDiagonal(int*, const std::string &, const std::string &, int*);

	inline bool QuickAlignCheck(int* algnPos, double *matchThresh, const std::string & query_seq, const std::string & cluster_seq){
		int scores[2] = { 1, 0 }; //match score is 1, mismatch score is 0, that we can count only matches
		//algnPos[1] = std::max(0, algnPos[1]);
		//algnPos[2] = std::min(algnPos[2], algnPos[3] - algnPos[0], algnPos[4] - algnPos[1]); //length of comparison will be limited by either the entire length of each sequence - start position, or by the desired length (input of algnPos[2])
		double numMatch = AlignDiagonal(algnPos, query_seq, cluster_seq, scores);
		//double numMatch = std::inner_product(query_seq.begin() + algnPos[0], query_seq.begin() + algnPos[0] + algnPos[2], cluster_seq.begin() + algnPos[1], 0, std::plus<int>(), compute_algn_score(scores[0], scores[1])); //count number of matches between sequences

		return numMatch > matchThresh[algnPos[2]];
	}

	int k_element(double *, int, int);

	inline bool Compare2DVec(const std::vector<double> & a, const std::vector<double> & b){
		return(a[0] > b[0]);
	}

	inline int AdjustToCorrectFrame(int queryStart, int germlineStart, int germline_reading_frame){
		int shift_start_pos = (germlineStart % 3 + 1) - germline_reading_frame,codonStart=queryStart;
		codonStart -= shift_start_pos;		
		while (codonStart < queryStart)
			codonStart += 3;			
		return codonStart;
	}

	bool CompareInt(int a, int b);

	void ShiftAlignPosition(int *, int, int);

	void CheckPeptideGapPositions(int *, int *, int *, int *, int ***);

	void FixGaps(std::string &, int **, int **, const int &, const int &);
}


#endif