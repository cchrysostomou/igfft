#include <iostream>
#include <algorithm> // min, max, min_element, max_element other stuff?
#include <fstream> //for ifstream
#include <cstring> //for strlen
#include <cstdio>
#include <vector>
#include <typeinfo> //for typeid
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string>
#include <sstream> //to use istringstream and getline
#include "StructVars.h"
#include <functional>
#include <cctype>
#include <locale>
#include <complex>
#include <set>
#include <map>
#include "DnaFunctions.h"
#include <stdexcept>      // std::out_of_range
using namespace std;

void dnafunctions::FixGaps(string & mod_seq, int **insertions, int **deletions, const int & numInsRow, const int & numDelRow){
	int insPoints = 0, delPoints = 0;

	for (int j = 0; j < numInsRow; j++)
	{
		mod_seq.insert(insertions[j][0] + insPoints + 1, insertions[j][1], NOT_NT); //insert command works by inserting to the LEFT of the  position provided. our insertions are defeined as "position right before gap"
		insPoints += insertions[j][1];
		for (int z = 0; z<numDelRow; z++)
		if (deletions[z][0]>insertions[j][0])
			deletions[z][0] += insPoints;
	}

	for (int i = 0; i < numDelRow; i++){
		mod_seq.erase(deletions[i][0] - delPoints, deletions[i][1]);
		delPoints += deletions[i][1];
	}

	return;
}

const map<std::string, char> dnafunctions::GenerateCodonMap(){
	map<string, char> temp_codon_map;

	set<char> unique_aa;
	string base;
	string codon(3, 'N');
	vector<char> all_bases;

	//MAP THE CODON TABLE OF NUCLEOTIDES TO AMINO ACID 
	//these are based on teh codon table and do not include ANY ambiguous bases 
	const map<string, char> default_codonmap = {
		//FIRST CODON: A
		//SECOND CODON: A				
		{ "AAA", 'K' },
		{ "AAC", 'N' },
		{ "AAG", 'K' },
		{ "AAT", 'N' },

		//SECOND CODON: C				
		{ "ACA", 'T' },
		{ "ACC", 'T' },
		{ "ACG", 'T' },
		{ "ACT", 'T' },

		//SECOND CODON: G				
		{ "AGA", 'R' },
		{ "AGC", 'S' },
		{ "AGG", 'R' },
		{ "AGT", 'S' },

		//SECOND CODON: T
		{ "ATA", 'I' },
		{ "ATC", 'I' },
		{ "ATG", 'M' },
		{ "ATT", 'I' },

		//FIRST CODON: C
		//SECOND CODON: A
		{ "CAA", 'Q' },
		{ "CAC", 'H' },
		{ "CAG", 'Q' },
		{ "CAT", 'H' },

		//SECOND CODON: C
		{ "CCA", 'P' },
		{ "CCC", 'P' },
		{ "CCG", 'P' },
		{ "CCT", 'P' },

		//SECOND CODON: G
		{ "CGA", 'R' },
		{ "CGC", 'R' },
		{ "CGG", 'R' },
		{ "CGT", 'R' },

		//SECOND CODON: T
		{ "CTA", 'L' },
		{ "CTC", 'L' },
		{ "CTG", 'L' },
		{ "CTT", 'L' },

		//FIRST CODON: G
		//SECOND CODON: A
		{ "GAA", 'E' },
		{ "GAC", 'D' },
		{ "GAG", 'E' },
		{ "GAT", 'D' },

		//SECOND CODON: C
		{ "GCA", 'A' },
		{ "GCC", 'A' },
		{ "GCG", 'A' },
		{ "GCT", 'A' },

		//SECOND CODON: G
		{ "GGA", 'G' },
		{ "GGC", 'G' },
		{ "GGG", 'G' },
		{ "GGT", 'G' },

		//SECOND CODON: T
		{ "GTA", 'V' },
		{ "GTC", 'V' },
		{ "GTG", 'V' },
		{ "GTT", 'V' },

		//FIRST CODON: T
		//SECOND CODON: A
		{ "TAA", '*' },
		{ "TAC", 'Y' },
		{ "TAG", '*' },
		{ "TAT", 'Y' },

		//SECOND CODON: C
		{ "TCA", 'S' },
		{ "TCC", 'S' },
		{ "TCG", 'S' },
		{ "TCT", 'S' },

		//SECOND CODON: G
		{ "TGA", '*' },
		{ "TGC", 'C' },
		{ "TGG", 'W' },
		{ "TGT", 'C' },

		//SECOND CODON: T
		{ "TTA", 'L' },
		{ "TTC", 'F' },
		{ "TTG", 'L' },
		{ "TTT", 'F' },

	};

	//maps ambigious DNA BASE calls to an array of potential bases 
	map<char, std::vector<char>> ambiguous_bases = {
		{ 'N', { 'A', 'C', 'G', 'T' } },
		{ 'X', { 'A', 'C', 'G', 'T' } },
		{ 'U', { 'T' } },
		{ 'K', { 'G', 'T' } },
		{ 'M', { 'A', 'C' } },
		{ 'R', { 'A', 'G' } },
		{ 'Y', { 'C', 'T' } },
		{ 'S', { 'C', 'G' } },
		{ 'W', { 'A', 'T' } },
		{ 'B', { 'C', 'G', 'T' } },
		{ 'V', { 'A', 'C', 'G' } },
		{ 'H', { 'A', 'C', 'T' } },
		{ 'D', { 'A', 'G', 'T' } },
		{ 'A', { 'A' } },
		{ 'C', { 'C' } },
		{ 'G', { 'G' } },
		{ 'T', { 'T' } },
	};

	for (map < char, vector<char>>::iterator new_base = ambiguous_bases.begin(); new_base != ambiguous_bases.end(); new_base++){
		//add ambiguous bases to all_bases call 
		all_bases.push_back(new_base->first);
	}

	//Generate list of all possible codons using bases listed in the vector of bases (all-bases)
	vector<string> codons;
	for (int i = 0; i < all_bases.size(); i++){
		codon[0] = all_bases[i];
		for (int j = 0; j < all_bases.size(); j++){
			codon[1] = all_bases[j];
			for (int k = 0; k < all_bases.size(); k++){
				codon[2] = all_bases[k];
				codons.push_back(codon);
			}
		}
	}

	char current_base;
	vector<vector<char>> potential_codons(3);
	potential_codons[0].push_back('a');
	vector<string> all_codon_search;
	//loop through all posssible codons determiend by fucntion above (possible codones = len(ambiguous_bases)^3)
	for (int i = 0; i < codons.size(); i++){
		//for each codon, we will look at how many unique amino acids we find 
		unique_aa.clear();
		codon = codons[i];
		for (int r = 0; r < 3; r++)
			potential_codons[r].clear();
		all_codon_search.clear();
		//loop through all th positions of the codon
		for (int pos = 0; pos < 3; pos++){
			current_base = codons[i][pos];
			//for each ambiguous base in the codon, loop through all of its potential bases based on ambigous bases map 
			for (vector<char>::iterator possible_base = ambiguous_bases[current_base].begin(); possible_base != ambiguous_bases[current_base].end(); possible_base++)
				potential_codons[pos].push_back(*possible_base);
		}

		for (int i1 = 0; i1 < potential_codons[0].size(); i1++){
			for (int j1 = 0; j1 < potential_codons[1].size(); j1++){
				for (int k1 = 0; k1 < potential_codons[2].size(); k1++){
					string c = "";
					c += potential_codons[0][i1];
					c += potential_codons[1][j1];
					c += potential_codons[2][k1];

					all_codon_search.push_back(c);
				}
			}
		}

		for (vector<string>::iterator it = all_codon_search.begin(); it != all_codon_search.end(); it++){
			//add all unqiue amino acids to teh set
			if (default_codonmap.find(*it) != default_codonmap.end()) //only will add amino acids found in codon map (so ignores ambiguous bases  from other positions, clunky code, but works)
				unique_aa.insert(default_codonmap.at(*it));
		}

		if (unique_aa.size() == 1)
			//this combination of bases (even ambiguous ones) results in a single unique amino acid. i.e AGN => SERINE
			//so add this codon to the codontable (any codon nto int he codon table will be reported as 'X')
			temp_codon_map[codons[i]] = *unique_aa.begin();
	}

	const map<string, char> actual_codon_map = temp_codon_map;

	return temp_codon_map;
}

//ntstart => the starting position of the nucloetide sequence (with regard tothe oroginal full length sequence)
//frame => the codon frame WITH REGARD TO THE ORIGINAL FULL LENGTH SEQUENCE (not nt_sequence) 
std::string dnafunctions::TranslateSeq(const string & nt_sequence,int frame,int ntstart){
	frame =  (frame - 1)-(ntstart%3);
	if (frame < 0)
		frame += 3;
	int aa_len = (nt_sequence.length() - frame);
	
	if (aa_len < 3 || ntstart == -1)
		return "";
		
	string nt = nt_sequence;
	string protein;
	
	string codon;	

	char aa;	
	int add_n = 3-(aa_len) % 3;
	if (add_n != 3)
		nt.append(add_n,'N');
	protein.resize(nt.length() / 3);

	map<string, char>::const_iterator it;
	int aa_pos = 0;
	int stop_len = nt.length() - 2;
	

	string substr;
	substr.resize(3);
	string::iterator chit;
	
	for (chit = nt.begin()+frame; chit<nt.end(); ++chit){
		substr[0] = *chit;
		chit++;
		substr[1] = *chit;
		chit++;
		substr[2] = *chit;
		it = dna_to_aa_map.find(substr);
		if (it == dna_to_aa_map.end()){			
			protein[aa_pos] = 'X';
		}
		else{
			protein[aa_pos] = it->second;
		}
		aa_pos += 1;
	}
				
	return protein.erase(protein.find_last_not_of('X')+1);
}

std::string dnafunctions::ReverseComplement(string & nt_sequence){
	string revComp = nt_sequence;
	char val;
	int a = 0;
	for (int i = nt_sequence.length() - 1; i >= 0; --i)
	{
		val = nt_sequence[i];
		switch (val)
		{
		case 'A':
		case 'a':
			revComp[a] = 'T';
			break;
		case 'T':
		case 't':
			revComp[a] = 'A';
			break;
		case 'G':
		case 'g':
			revComp[a] = 'C';
			break;
		case 'C':
		case 'c':
			revComp[a] = 'G';
			break;
		default:
			nt_sequence[i] = NOT_NT;// char(0);
			revComp[a] = NOT_NT;// char(0);
			break;
		}
		a = a + 1;
	}
	return revComp;
}

double dnafunctions::AlignDiagonal(int *info, const std::string & sequence, const std::string & germline, int *algnScores){
	int lenCompare = (info[0] + info[2] - 1 > sequence.length() - 1 || info[1] + info[2] - 1 > germline.length() - 1) ? std::min(sequence.length() - info[0], germline.length() - info[1]) : info[2];// min(sequence.length() - info[0] , germline.length() - (info[1] + info[2]) - 1);

	int matches[3] = { 0 }; //matches[1] = number match bases, matches[0] = number mismatching bases  nMatch = 0, nMismatch = 0; maches[2] => no penalty or match, represnts x match with x
	int compare;

	for (int i = 0; i < lenCompare; i++){
		compare = (sequence[i + info[0]] == NOT_NT || germline[i + info[1]] == NOT_NT) ? 2 : (sequence[i + info[0]] == germline[i + info[1]]);
		matches[compare]++;
	}

	double score = matches[1] * algnScores[0] - matches[0] * (abs(algnScores[1]));
	info[2] = matches[1] + matches[0];
	return score;
}

int dnafunctions::k_element(double *v, int numfind, int length){
	std::sort(v, v + length);
	int numHit = 1;
	double max_hit;
	int i = length - 1;
	bool loop;
	do{
		i--;
		max_hit = v[i + 1];
		if (v[i] != max_hit){
			numHit++;
		}
		loop = (i < 0 || numHit > numfind);
	} while (!loop);
	return i + 1;
}

void dnafunctions::ShiftAlignPosition(int *shiftedAlgnPos, int seqStart, int germEnd){
	int shift;
	//int algnP2 = shiftedAlgnPos[1];
	//int algnP4 = shiftedAlgnPos[4];
	shiftedAlgnPos[1] -= 1; //move the index of sequence/character position in string down by one to match index[0] = position one
	shiftedAlgnPos[2] -= 1;
	shiftedAlgnPos[3] -= 1;
	shiftedAlgnPos[4] -= 1;

	shift = seqStart > shiftedAlgnPos[1] ? seqStart - shiftedAlgnPos[1] : 0; //if the starting position of the sequence is at a position higher than where the alignment to the querys equence begins, then shift the alignment to start where they strting squence is
	shiftedAlgnPos[1] += shift;
	shiftedAlgnPos[3] += shift;

	shift = germEnd < shiftedAlgnPos[4] ? germEnd - shiftedAlgnPos[4] : 0; //if the alignment position is a position greater than the end of the germline, then readjust position
	shiftedAlgnPos[2] += shift;
	shiftedAlgnPos[4] += shift;

	shiftedAlgnPos[5] = abs(shiftedAlgnPos[2] - shiftedAlgnPos[1] + 1);

	return;
}

/*%once we perform the FFT alignment to identify the best "diagonal" or
%maximu region of alignment, we cannot be sure if there might be a few
%gaps at the beginning or end of the sequence.
%so this functions, uses the region of multiple alignment to try to
%identify bases at the top and bottom to align to the germline to see if a
%gap ocurred upstream/downstream*/

//algnPos - > information about the alignment between the sequence query and germline
//peptideInfo -> position 0/first index = length of peptide, position1/second index = length of a gap. length of sequence allowed around gap
//seqInfo -> seq start/seq end => position 0 and 1 respectively
//germInfo -> germ start/germ end => position 0 and 1 respectively

void dnafunctions::CheckPeptideGapPositions(int *algnPos, int *peptideInfo, int *seqInfo, int *germInfo, int ***algnPepResults){
	//initialize algnPepResults
	// each row in algnPepResuts are as such: begP, endP, begP_ds50, endP_us50, midP
	//each value is a 2x2 array of [0,0;0,0]
	for (int i = 0; i < 5; i++){
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){
				algnPepResults[i][j][k] = 0;
			}
		}
	}

	//algnPepResults[5] = 0; //algnPepResults[5] => corrrespondes to algnLen

	int algnLen = algnPos[2] - algnPos[1];
	int check2dist = 50;
	int startP = algnPos[1] - peptideInfo[1]; //algnPos[1] - gapSeq
	int stopSequence;

	if (startP < seqInfo[0]) //startP<seqStart
		algnPepResults[0][0][0] = seqInfo[0];//begP(1,1) = seqStart
	else
		algnPepResults[0][0][0] = startP; //begP(1,1) = startP

	algnPepResults[0][1][0] = algnPos[3];//begP(2,1) = algnPos[3]
	algnPepResults[0][0][1] = algnPepResults[0][0][0] + peptideInfo[0] + 2 * peptideInfo[1] - 1; //begP(1,2) = begP(1,1)+peptideSize+2*gapSeq-1
	algnPepResults[0][1][1] = algnPepResults[0][1][0] + peptideInfo[0] - 1; //begP(2,2) = begP(2,1)+peptideSize-1

	if (algnLen>check2dist){
		stopSequence = algnPos[2] + peptideInfo[1]; //stopSequence = algnPos[2]+gapseq
		if (stopSequence > seqInfo[1]) //if stopsequence>seqEnd
			algnPepResults[1][0][1] = seqInfo[1]; //endP(1,2) = seqEnd
		else
			algnPepResults[1][0][1] = stopSequence; //endP(1,2) = stopsequence
		if (algnPos[4] > germInfo[1]) // if algnPos[4] > germEnd
			algnPepResults[1][1][1] = germInfo[1]; //endP(2,2) = germEnd
		else
			algnPepResults[1][1][1] = algnPos[4]; //endP(2,2) = algnPos[4]
		algnPepResults[1][1][0] = algnPepResults[1][1][1] - peptideInfo[0] + 1; //endP(2,1) = endP(2,2)-peptideize+1
		algnPepResults[1][0][0] = algnPepResults[1][0][1] - 2 * peptideInfo[1] - peptideInfo[0] + 1; //endP(1,1) = endP(1,2)-2*gapSeq-peptideSize+1
	}

	double middleNTPos;
	//if the sequence alignment is longer than 100 bp , then we will report the
	if (algnLen > check2dist * 2){
		middleNTPos = round((algnLen / 2.0) - (peptideInfo[0] / 2.0));  //middleNTPos = round(algnLen/2 - peptidesize/2)
		algnPepResults[4][0][0] = algnPepResults[0][0][0] + middleNTPos - 1; //midP(1,1) = begP(1,1)+middleNTPos-1
		algnPepResults[4][0][1] = algnPepResults[0][0][1] + middleNTPos - 1;//midP(1,2) = begP(1,2)+middleNTPos-1
		algnPepResults[4][1][0] = algnPepResults[0][1][0] + middleNTPos - 1;//midP(2,1) = begP(2,1)+middleNTPos-1
		algnPepResults[4][1][1] = algnPepResults[0][1][1] + middleNTPos - 1;//midP(2,2) = begP(2,2)+middleNTPos-1
	}

	//%look at positions about 50 bases (defined by check2Dist) downstream of
	//%alignment start
	int temp1, temp2, temp3, temp4;
	if (algnLen > check2dist * 3){
		temp1 = algnPepResults[0][0][1] + check2dist - 1; //temp1 = begP(1,2) +check2dist-1;
		temp2 = algnPepResults[0][1][1] + check2dist - 1; //temp2 = begP(2,2) +check2dist-1;

		if (temp1 <= seqInfo[1] && temp2 <= germInfo[1]){ // %first make sure that the nucleotide positions of the ds50 sequence is actualy within the alignment range
			algnPepResults[2][0][1] = temp1;//begP_ds50 = temp1  // %temp1 was just the end of the sequence position for this germline
			algnPepResults[2][0][0] = algnPepResults[0][0][0] + check2dist - 1; //begP_ds50(1,1) = begP(1,1)+check2dist-1
			algnPepResults[2][1][0] = algnPepResults[0][1][0] + check2dist - 1; //begP_ds50(2,1) = begP(2,1)+check2dist-1
			algnPepResults[2][1][1] = temp2;// %temp2 was just the end of the germline position for this peptide
		}

		temp3 = algnPepResults[1][0][0] - check2dist + 1; //temp3 = endP(1,1)-check2dist+1
		temp4 = algnPepResults[1][1][0] - check2dist + 1; //temp4 = endP(2,1) -chekc2dist+1

		if (temp3 >= seqInfo[0] && temp4 >= germInfo[0]){ //if temp3>=seqstart && temp4>=germStart
			algnPepResults[3][1][0] = temp4; //enP_us50(2,1) = temp4
			algnPepResults[3][1][1] = algnPepResults[1][1][1] - check2dist + 1; //endP_us50(2,2) = endP(2,2)-check2dist+1
			algnPepResults[3][0][0] = temp3;
			algnPepResults[3][0][1] = algnPepResults[1][0][1] - check2dist + 1;//endP_us50(1,2) = endP(1,2) -chekc2dist+1
		}
	}

	return;
}