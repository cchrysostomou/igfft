#include <iostream>
#include <iomanip>
#include "string.h"
#include "DnaFunctions.h"
#include "StructVars.h"
using namespace std;

int main(){
	int swAlignCoordinates[5];
	string seqA = "GAGGTGCAGCTGGTGGAGTCCGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGGTTCACCTTCAGTGGCCCCGCCATGCACTGGGTCCGCCAGGCTCCCGGGAAGGGGCTGGAGTGGGTTGGCCGTACTAGAAGCAAAGCTAACAGTTACGCCACAGCATACGCCGCGTCGGTGAAAGGCAGGTTCACCATCTCCAGAGATGATTCAAAGAACACGCCGTATCTGCAAATGAACAGCCTGAAAACCGAGGACACGGCCGTGTATTACTGTACTAGACA"; //query
	string seqB = "GAGGTGCAGCTGGTGCAGTCTGGAGCAGAGGTGAAAAAGCCCGGGGAGTCTCTGAGGATCTCCTGTAAGGGTTCTGGATACAGCTTTACCAGCTACTGGATCGGCTGGGTGCGCCAGATGCCCGGGAAAGGCCTGGAGTGGATGGGGAGCATCGATCCTGGTGACTCTGATACCAGCTACAGCCCGTCCTTCCAAGGCCAGGTCACCATCTCAGCCGACAAGTCCATCAGCACCGCCTACCTGCAGTGGAGCAGCCTGAAGGCCTCGGACACCGCCCCGTATTACTGTGCGCGACA"; //germline

	swAlignCoordinates[0] = 0;
	swAlignCoordinates[1] = seqB.length();
	swAlignCoordinates[2] = 0;// best_cluster_so_far->alignment_to_sequence[1] - best_cluster_so_far->alignment_to_sequence[3]; //diagonal defined as sequence - germline
	swAlignCoordinates[3] = 0;
	swAlignCoordinates[4] = 6;

	structvars::SWAlignment aligner(seqA.length(), seqB.length(), 10, 5, 4, 50, 10);
	dnafunctions::SWAlignDiag(seqA, seqB, swAlignCoordinates, aligner);
	dnafunctions::TraceBack(seqA, seqB, aligner);

	cout << "stop" << endl;
	return 0;
}