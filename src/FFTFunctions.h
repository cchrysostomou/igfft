#ifndef _FFT_FUNCTIONS_H
#define _FFT_FUNCTIONS_H
//#include <ctime>
#include "StructVars.h"
#include "fftw3.h"
#include <map>
//#include <sys/time.h>
#include <time.h>

namespace ffthelper{
	//bool IsDivisible(int); //check if value is divisible by predetermined prime numbers
	bool IsDivisible(int, std::vector<int>, std::vector<int> &);
	void makeFakeSeqInt(int n, complex_t *);
	//void AnalyzeFFTTime(int, int, double*);
	double AnalyzeFFTTime(int); //measure the number of FLOPS for an FFTW plan given the length of n (int)

	//the fft function is optimized for certain lengths.  In Theory it should scalle as NlogN (where log is base 2)
	// occasionally we will hit an alignment size that is not the "best" choice for the fftw program (i.e. 791 apparently is a bad length to perform fft)
	//so this function will make sure the length we use is indeed NLogN based on the "base 2" weighting system
	int FindIdealAlgnSize(int);
	
	void MapCrossCorrelationToShift(int **&model, int sizeModel, const std::vector<int> & seq1, const std::vector<int> & seq2);

	//////////INLINE FUNCTIONS/////////////////
	inline void PadSeqZerosEnd(complex_t *dnaInt, complex_t *paddedSeq, int seqLength, int paddedLength){
		//(add zeroes to the end)
		memcpy(paddedSeq, dnaInt, (seqLength)*sizeof(complex_t)); //copy contents of one to the other
		for (int i = seqLength; i < paddedLength; i++){
			paddedSeq[i][0] = 0;
			paddedSeq[i][1] = 0;
		}
	}

	inline void PadSeqZerosFront(complex_t *dnaInt, complex_t *paddedSeq, int seqLength, int paddedLength){
		//copy elements to the end of the pointer (Add zeros to the front)
		int byteToMove = paddedLength - seqLength;
		complex_t *tmp = (complex_t*)malloc((seqLength)*sizeof(complex_t)); //make temp array
		memcpy(tmp, dnaInt, (seqLength)*sizeof(complex_t)); //copy contents from original to temp
		memmove(paddedSeq + (byteToMove), tmp, (seqLength)*sizeof(complex_t)); //move the new pointer already padded
		for (int i = 0; i < byteToMove; i++){
			paddedSeq[i][0] = 0;
			paddedSeq[i][1] = 0;
		}
		free(tmp);
	}


	//void convertSeq(complex_t*, std::string &); //function: convert a dna sequence into a set of complex integers
	//function: convert a dna sequence into a set of complex integers
	inline void convertSeq(complex_t *dnaInt, std::string &dnaSeq)
	{
		//map<char, int> letConv = { { 'A', 1 }, { 'T', 2 }, { 'G', 3 }, { 'C', 4 }, { 'a', 1 }, { 't', 2 }, { 'g', 3 }, { 'c', 4 } };
		int seqSize = dnaSeq.length();
		//method - we will use this to change the different ways in which we can digitize the dna sequence
		//result = (complex_t *)calloc(seqSize, sizeof(complex_t));
		for (int i = 0; i < seqSize; ++i){//for letter in the sequence create the "integer" representation
			switch (dnaSeq[i])
			{
			default:
				dnaInt[i][0] = 0;
				dnaInt[i][1] = 0;
				break;
			case 'A':
			case 'a':
				dnaInt[i][0] = 0;
				dnaInt[i][1] = 1;
				break;
			case 'T':
			case 't':
				dnaInt[i][0] = 0;
				dnaInt[i][1] = -1;
				break;
			case 'G':
			case 'g':
				dnaInt[i][0] = -1;
				dnaInt[i][1] = 0;
				break;
			case 'C':
			case 'c':
				dnaInt[i][0] = 1;
				dnaInt[i][1] = 0;
				break;
			}
		}
	}

	inline void CrossCorrelation(complex_t *cross_vector, complex_t *germline, complex_t *seq, int n){
		for (int i = 0; i < n; i++){
			//alignment.crossCorrVector[i][0] = seqSignal[i][0]*db_signal[i][0];
			//alignment.crossCorrVector[i][1] = seqSignal[1][0]*db_signal[i][0];
			cross_vector[i][0] = (seq[i][0] * germline[i][0]) - (seq[i][1] * germline[i][1]);
			cross_vector[i][1] = (seq[i][1] * germline[i][0]) + (seq[i][0] * germline[i][1]);
		}
	}
}

#endif