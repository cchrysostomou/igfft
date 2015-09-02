//THIS PROGRAM RUNS A TEST TO PERFORM FFT FOR DNA SEQUENCE ALIGNMENTS!!

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "math.h"
#include <time.h>
#include <algorithm>
#include <string>
#include <cmath>
#include <time.h>

#include <complex.h> //always include complex.h before including fftw3.  thilw will use C-version of complex numbers (not c++???).  I think its more effiicient because data shows twice improvement of speed with using this 
#include "fftw3.h" //this is the only library needed for using FFTW 


using namespace std;

//function: convert a dna sequence into a set of complex integers
fftw_complex *convertSeq(string dnaSeq, int defSize, int method)
{
	//method - we will use this to change the different ways in which we can convert seqeuence to series of cmplex numbers
	int defaultSize = defSize;
	
	fftw_complex *dnaInt = (fftw_complex *)calloc(defaultSize, sizeof(*dnaInt));

	return dnaInt;
}



//function: convert DNA sequences into integers and fourier transform
void convertInteger(vector<seqDat> &sequence, int debug)
{
	int defSize = 512;
	string dnaSeq;
	fftw_complex *x = (fftw_complex*)calloc(defSize, sizeof(*x));
	fftw_complex *xf = (fftw_complex*)calloc(defSize, sizeof(*xf));
	fftw_complex val;

	for (int i = 0; i<sequence.size(); ++i)
	{
		xf = (fftw_complex*)calloc(defSize, sizeof(*xf));
		dnaSeq = sequence[i].seq;
		sequence[i].seqInt = convertSeq(dnaSeq, defSize, 1); //create the integer value of the sequence
		x = sequence[i].seqInt;
		fftw_dft(xf, defSize, x); //calculate FFT of DNA sequence
		sequence[i].fftInt = xf; //store result in struct	
	}
}

//Haitham: this is the core of the function, it performs the cross correlation for each alignment and then the inverse FFT
//this is what takes 12 minutes to run


//Function: uses FFT/Cross correlation to identify the best match of each sequence against the possible germline sequences
void findBestHit(vector<seqDat> sequence, vector<seqDat> germline, vector<seqDat> germRev, int n)
{
	fftw_complex *seq;
	fftw_complex *gF;
	fftw_complex *gR;
	fftw_complex *cross = (fftw_complex*)calloc(n, sizeof(*cross));
	fftw_complex *crossR = (fftw_complex*)calloc(n, sizeof(*crossR));
	fftw_complex *x = (fftw_complex*)calloc(n, sizeof(*x));
	fftw_complex *xr = (fftw_complex*)calloc(n, sizeof(*xr));

	double fwdA[n];
	double revA[n];
	double *maxHit = 0;
	double *maxHitR = 0;
	int bestGene = 0;

	fftw_plan p = fftw_plan_dft_1d(n, cross, x, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan pr = fftw_plan_dft_1d(n, crossR, xr, FFTW_BACKWARD, FFTW_ESTIMATE);


	vector<results> summary(sequence.size());
	results temp;

	double t1 = 0, t2 = 0, t3 = 0, t4 = 0, d;
	reset_timer();

	for (int i = 0; i<sequence.size(); ++i)
	{
		if (i % 5000 == 0){
			cout << i << "----" << n << "----" << sequence.size() << "----" << germline.size() << "----" << t1 << "----" << t2 << "----" << t3 << "----" << t4 << endl;
			t1 = 0;
			t2 = 0;
			t3 = 0;

		}

		temp.maxScore = 0;
		temp.dir = 'f';

		//loop through every sequence where we are trying to associate a V gene with it
		//retrieve the stored fast fourier for this sequence
		seq = sequence[i].fftInt;

		//loop through every possible germline gene (we are aligning each sequence to each possible germline gene)
		for (int j = 0; j<germline.size(); ++j)
		{
			//retrieve the stored fast fourier for the forward diretion of this germline Gene
			gF = germline[j].fftInt;

			//retrieve the stored fast fourier for the reverse direction of this germline Gene (remember double stranded DNA, two possible signals)
			gR = germRev[j].fftInt;

			//               d=get_time();
			for (int k = 0; k<n; ++k){
				//perform cross correlation of sequence to forward and reverse DNA sequences
				cross[k] = seq[k] * conj(gF[k]); //sliding dot product = direct multiplication in FFT space
				crossR[k] = seq[k] * conj(gR[k]);

			}

			//           t1 +=get_time() - d;
			//             d= get_time();

			//calculate iFFT of cross correlations
			fftw_execute(p);
			fftw_execute(pr);
			//fftw_dft(x,n,cross,1); //calculate FFT of DNA sequence
			//fftw_dft(xr,n,crossR,1); //calculate FFT of DNA sequence

			//         t2 +=get_time() - d;
			//       d= get_time();


			//extract only the real values
			for (int k = 0; k<n; ++k){
				//			x[k]/=n;
				//			xr[k]/=n;
				fwdA[k] = creal(x[k]);
				revA[k] = creal(xr[k]);

			}

			//     t3 +=get_time() - d;
			//   d= get_time();

			//find the maximum alignment score in each cross-correlation
			maxHit = max_element(fwdA, fwdA + n);
			maxHitR = max_element(revA, revA + n);


			//choose the maximum of the two possible alignments 
			//sequence aligned to the forward germline DNA sequence
			//OR sequence aligned to the reverse of the germline DNA sequence, choose the beset
			if (*maxHit >= *maxHitR)
			{
				temp.dir = 'f';
				if (*maxHit>temp.maxScore){
					temp.vgene = germline[j].seqHeader; //store the name of the best germline hit
					temp.maxScore = *maxHit;
				}

			}
			else
			{
				temp.dir = 'r';
				if (*maxHitR>temp.maxScore){
					temp.vgene = germline[j].seqHeader;
					temp.maxScore = *maxHitR;
				}
			}

			//t4 +=get_time() - d;
			// d= get_time();							

		}

		summary[i] = temp; //summarize the results for each sequence
		summary[i].seq = sequence[i].seq;
		summary[i].seqHeader = sequence[i].seqHeader;
	}

	fftw_destroy_plan(p);
	fftw_destroy_plan(pr);

	cout << sequence.size() << "----" << n << "----" << sequence.size() << "----" << germline.size() << "----" << t1 << "----" << t2 << "----" << t3 << "----" << t4 << endl;
	//cout<<"TOTAL TIME:"<<get_time()<<endl;


}


