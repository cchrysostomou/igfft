#define _CRT_SECURE_NO_DEPRECATE

#ifndef FFTALIGN_H
#define FFTALIGN_H

#include <iomanip>
#include <iostream>
#include <assert.h>
#include <stdio.h>

#include <fstream>
#include "StructVars.h"
#include "FFTFunctions.h"

class FFTAlign{
	friend class QueryAlignment;
	friend class GermlineCluster;
public:
	FFTAlign(int);
	FFTAlign(int, structvars::FFTSettings);
	FFTAlign(int, const std::vector<int> &, const std::vector<int> &);
	FFTAlign(int, const std::vector<int> &, const std::vector<int> &, structvars::FFTSettings);
	
	FFTAlign(const FFTAlign &); //copy constructor
	~FFTAlign();

	void UpdateOverlapModel(const std::vector<int> &, const std::vector<int> &);
	double ReturnBinomialError(int, int); //first int is the row number, second int is the column number (alignment length in this case)
	double ReturnAlignmentScoreAtPosition(int);
	int GetAlignmentLength();

	void ForwardFFT(complex_t*, complex_t*);//coverts a sequence to forward fft
	int AlignSequences(complex_t *, complex_t *); //actually performs an fft alignment
	void IdentifyAdditionalAlignments(int, int col_num = 5);
	bool FindGappedPeaks(int *, double);

private:
	void CreateFFTPlan(); //when FFTalign is created, this function is called to create the proper variables, alignment length, and fft plans for performing transfrom
	void ModelBinomialAlignmentError(); //when FFTalign is created, this function will create a model of the expected cross correlation score given a random alignment length

	structvars::FFTSettings fftparams; //contains the following parameters: sensitivity, gapopen, gapextend
	int alignmentLength;
	int numGaps, bestDiag, minDiag, maxDiag;
	double gapPenalty, gappedScore, bestDiagScore, totalScore;
	complex_t *inputVector, *fourierVector, *crossCorrVector, *alignmentVector;
	int **alignmentModel;
	double **expectedErrorModel;
	fftw_plan forwardT, reverseT;
};

inline void FFTAlign::ForwardFFT(complex_t *complex_seq, complex_t *fourier_seq){//coverts a sequence to forward fft
	memcpy(inputVector, complex_seq, sizeof(complex_t)*alignmentLength); //copy contents to vector used for fftw
	fftw_execute(forwardT);//perform forward transformation
	memcpy(fourier_seq, fourierVector, sizeof(complex_t)*alignmentLength);//copy results back into fourier signal
	
	//std::ofstream testfile("fftval.txt");
	/*for (int i = 0; i < alignmentLength; i++)
		testfile << inputVector[i][0] << "\t" << inputVector[i][1] << "\n";
	testfile.close();*/
}

#endif