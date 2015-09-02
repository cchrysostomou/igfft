#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include "fftw3.h"
#include <vector>
#include <time.h>

typedef fftw_complex complex_t;

using namespace std;

void makeFakeSeqInt(int n, complex_t *x){
	srand(time(NULL));
	for (int i = 0; i<n; i++){
		int val = rand() % 4;
		switch (val){
		case 0:	x[i][0] = 0;
			x[i][1] = 1;
			break;
		case 1:
			x[i][0] = 0;
			x[i][1] = -1;
			break;
		case 2:
			x[i][0] = 1;
			x[i][0] = 0;
			break;
		case 3:
			x[i][0] = -1;
			x[i][1] = 0;
			break;
		}
	}

}

//22530000
// 4140000
// 3900000
// 5010000

int main(){
	int numRepeats = 19;
	//	int n = 1024;
	int numSeqs = 10000;

	vector<complex_t*> testVals;


	int numData = 1024;
	int fftw_init_threads(void);
	int nthreads = 4;
	//fftw_plan_with_nthreads(nthreads);	
	//fftw_plan test = fftw_plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_PATIENT);

	int rank = 1;
	int n[1] = { numData };
	int howmany = numRepeats * 2;
	const int *inembed = NULL;
	const int *onembed = NULL;
	int istride = 1; //location of the next element in a transform
	int ostride = 1; //location of the next element in a transform
	int idist = numData; //location of the next transform
	int odist = numData; //location of the next transform
	int sign = FFTW_FORWARD;
	unsigned flags = FFTW_PATIENT;



	//	fftw_plan 1dplan = fftw_plan_dft_1d(n,
	//	a = 0;



	vector< vector<complex_t*> > seqs(numSeqs);
	for (int i = 0; i<numSeqs; i++){
		complex_t* temp = (complex_t*)malloc(numData*sizeof(complex_t));
		complex_t* temp2 = (complex_t*)malloc(numData*sizeof(complex_t));
	
		seqs[i].push_back(temp);
		seqs[i].push_back(temp2);
		makeFakeSeqInt(numData, seqs[i][0]);
		makeFakeSeqInt(numData, seqs[i][1]);
	}

	for (int i = 0; i<numRepeats; i++){
		complex_t* temp3 = (complex_t*)malloc(numData*sizeof(complex_t));
		testVals.push_back(temp3);
		makeFakeSeqInt(numData, testVals[i]);
	}


	complex_t *x1 = (complex_t*)fftw_malloc(numData*sizeof(complex_t));
	complex_t *xmany = (complex_t*)fftw_malloc(2 * numRepeats*numData*sizeof(complex_t));
	complex_t *y1 = (complex_t*)fftw_malloc(numData*sizeof(complex_t));
	complex_t *ymany = (complex_t*)fftw_malloc(2 * numRepeats*numData*sizeof(complex_t));

	fftw_plan plan = fftw_plan_many_dft(rank, n, howmany,
		xmany, inembed,
		istride, idist,
		ymany, onembed,
		ostride, odist,
		sign, FFTW_PATIENT);

	complex_t aa;
	aa[0] = 1;
	aa[1] = 0;
	complex_t* b = (complex_t*)malloc(10*sizeof(complex_t));
	for (int l = 1; l < numSeqs; l++){
		int a = 0;
		for (int i = 1; i < numData; i++){
			for (int j = 1; j < numRepeats; j++){
				for (int k = 1; k < 2; k++){
					xmany[a][0] = seqs[l][k][i][0] * testVals[j][i][0];
					xmany[a][1] = seqs[l][k][i][1] * testVals[j][i][1];
					a++;
				}
			}
		}
		fftw_execute(plan);
	}
	


	//fftw_plan plan = fftw_plan_dft_1d(30,x,y,FFTW_FORWARD,FFTW_PATIENT);
	/*	for (int i = 0;i<30;i++){
	x[i][0] = i;
	x[i][1] = 0;
	}*/
	/*
	clock_t start = clock()
	for (int i = 0;i<100000;i++){
	for (int j = 0;j<19;j++){
	makeFakeSeqInt(n,x);
	fftw_execute(test);
	}
	}

	clock_t end = clock();
	*/

	fftw_execute(plan);
	fftw_destroy_plan(plan);

	//fftw_cleanup_threads();

	/*	cout << x[12][0] << endl;
	cout << x[12][1] << endl;
	cout << y[11][0] << endl;
	cout << y[11][1] << endl;*/
	cout << "done" << endl;
	//cout<<end-start<<endl;
	cin.get();
	return 0;
}

// g++ testmultithread.cpp -lfftw3_threads -lfftw3 -lm -lpthread

