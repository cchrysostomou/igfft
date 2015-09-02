#include "StructVars.h"
#include "FFTFunctions.h"

#include <complex>
#include "fftw3.h"

using namespace std;
using namespace ffthelper;

//function: convert a dna sequence into a set of complex integers
complex_t *ffthelper::convertSeq(string dnaSeq, int defSize)
{
	//method - we will use this to change the different ways in which we can digitize the dna sequence
	int defaultSize = defSize;
	complex_t *dnaInt = (complex_t *)calloc(defaultSize, sizeof(*dnaInt));
	//complex_t dnaInt[defaultSize];
	double *each = (double*)calloc(50, sizeof(double*));
	complex_t I = { 0, 1 };
	
	char val;
	for (int i = 0; i<dnaSeq.length(); ++i){//for letter in the sequence create the "integer" representation
		val = dnaSeq[i];
		switch (val)
		{
		case 'A':
		case 'a':			
			dnaInt[i] = 1;
			break;
		case 'T':
		case 't':
			dnaInt[i] = -1;
			break;
		case 'G':
		case 'g':
			dnaInt[i] = I;
			break;
		case 'C':
		case 'c':
			dnaInt[i] = -I;
			break;
		default:
			dnaInt[i] = 0;
			break;
		}
	}	
	//cout << dnaInt[1] <<" "<< dnaInt[2] <<" "<< dnaInt[3] <<" "<< dnaInt[4] <<" "<< dnaInt[5]<< endl;
	for (int i = dnaSeq.length(); i<defaultSize; ++i)
		dnaInt[i] = 0;

	return dnaInt;
}


/*
//function: convert DNA sequences into integers and fourier transform
void ffthelper::convertInteger(vector<seqDat> &sequence, int debug)
{
	int defSize = 512;
	string dnaSeq;
	complex_t *x = (complex_t*)calloc(defSize, sizeof(*x));
	complex_t *xf = (complex_t*)calloc(defSize, sizeof(*xf));
	complex_t val;

	for (int i = 0; i<sequence.size(); ++i)
	{
		xf = (complex_t*)calloc(defSize, sizeof(*xf));
		dnaSeq = sequence[i].seq;
		sequence[i].seqInt = convertSeq(dnaSeq, defSize, 1); //create the integer value of the sequence
		x = sequence[i].seqInt;
		fftw_dft(xf, defSize, x); //calculate FFT of DNA sequence
		sequence[i].fftInt = xf; //store result in struct	
	}

}*/

//Initializing/making an FFTW plan for the alignment.  
void ffthelper::CreateFFTPlan(structvars::GermlineAlignmentSettings & algn){
	int n = algn.alignmentLength;
	algn.fT = fftw_plan_dft_1d(n, algn.inputVector, algn.fourVector,  FFTW_FORWARD, FFTW_ESTIMATE);
	algn.rT = fftw_plan_dft_1d(n, algn.crossCorrVector, algn.alignmentVector, FFTW_FORWARD, FFTW_ESTIMATE);	
}
	