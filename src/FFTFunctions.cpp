#include "FFTFunctions.h"

using namespace std;
using namespace ffthelper;

void ffthelper::makeFakeSeqInt(int n, complex_t *x){
	srand(time(NULL));
	for (int i = 0; i < n; i++){
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

//old and slow function
/*void ffthelper::AnalyzeFFTTime(int n, int numLoops, double *results){
	complex_t *x = (complex_t*)malloc(n*sizeof(complex_t));
	complex_t *y = (complex_t*)malloc(n*sizeof(complex_t));
	complex_t *xFFT = (complex_t*)malloc(n*sizeof(complex_t));
	fftw_plan setup = fftw_plan_dft_1d(n, x, xFFT, FFTW_FORWARD, FFTW_PATIENT);
	fftw_plan runFFT = fftw_plan_dft_1d(n, xFFT, y, FFTW_BACKWARD, FFTW_PATIENT);
	makeFakeSeqInt(n, x); //make a fake DNA sequence
	fftw_execute(setup); //convert to fourier singal
	clock_t start, end, g, time;

	double mean = 0, mean_sq = 0;

	start = clock();
	for (int i = 0; i < numLoops; i++)
		fftw_execute(runFFT); //test inverse
	end = clock();
	time = end - start;

	results[0] = time / (1.0*numLoops);
	results[1] = time*1.0;

	fftw_destroy_plan(runFFT);
	std::free(x);
	std::free(y);
	std::free(xFFT);
}*/

/*
int ffthelper::FindIdealAlgnSize(int target_len){
int numTest = 100000 / target_len;
int base = int(log2(target_len));
int lowVal = pow(2, base);
int highVal = pow(2, base + 1); //This is the worst case scenario we will ever consider (padding the data to the next highest power of 2)
double timeLow[2], timeHigh[2], timeVal[2];

//First lets figure out the worst case scenario . that is we have to pad n to the value highVal
AnalyzeFFTTime(highVal, numTest, timeHigh);
AnalyzeFFTTime(lowVal, numTest, timeLow);
double minTime = timeHigh[0];
double bestHit = highVal;
for (int i = target_len; i < highVal; i++){
if (IsDivisible(i)){ //if the value is divisible by the prime number
AnalyzeFFTTime(i, numTest, timeVal); //if its a good number, then lets analyze its speed, and check its indeed faster than being rounded to the next base 2 number
if (timeVal[0] <= minTime && i < bestHit){
minTime = timeVal[0];
bestHit = i;
}
}
}
return bestHit;
}
*/
/*
bool ffthelper::IsDivisible(int Val){

vector<int> ListOfAllowedPrimes = { 2, 3, 5, 7 };//analysis of fftw showed that the transform of only numbers divisible by these prime numbers result in proper computation time, The best possible computation time are those numbers divisible by 2 (but 3,5, and 7 are good too and sometimes 11 is good)
int n = Val;
for (int i = 0; i < ListOfAllowedPrimes.size(); i++){
while (n%ListOfAllowedPrimes[i] == 0)
n = n / ListOfAllowedPrimes[i];
}
if (n == 1)
return true;
else
return false;
}*/

//returns whether a function is divisible by the supplied vector (primes), and in addition, returns how each prime number factors into n (the powers vector) i.e. primes(1)^x*primes(2)^y = n where x and y = powers(1) and powers(2)
bool ffthelper::IsDivisible(int n, std::vector<int> primes, std::vector<int> & powers){

	powers.resize(primes.size());

	for (int i = 0; i < primes.size(); i++){
		powers[i] = 0;
		while (n%primes[i] == 0){
			n = n / primes[i];
			powers[i] += 1;
		}
	}
	if (n == 1)
		return true;
	else
		return false;
}


double ffthelper::AnalyzeFFTTime(int n){
	complex_t *testVec, *outputVec;
	fftw_plan test;
	double total;
	double add[1], mul[1], fma[1];
	testVec = (complex_t*)malloc(n * sizeof(complex_t));
	outputVec = (complex_t*)malloc(n * sizeof(complex_t));
	test = fftw_plan_dft_1d(n, testVec, outputVec, FFTW_FORWARD, FFTW_PATIENT);
	fftw_flops(test, add, mul, fma); //MEASURE THE NUMBER OF FLOPS FOR A PLAN OF SIZE N	
	total = *add + *mul + *fma * 2;//MEASURE THE NUMBER OF FLOPS FOR A PLAN OF SIZE N	

	fftw_destroy_plan(test);

	std::free(testVec);
	std::free(outputVec);

	return total;
}

//Find a length, N>=target_len, that minimized the number of computations required to perform a fourier transform
//this should optimize speed of transforms/gapless alignments
int ffthelper::FindIdealAlgnSize(int target_len){
	vector<int> ListOfAllowedPrimes = { 2, 3, 5, 7 };//analysis of fftw showed that the transform of only numbers divisible by these prime numbers result in proper computation time, The best possible computation time are those numbers divisible by 2 (but 3,5, and 7 are good too and sometimes 11 is good)
	int base = int(log2(target_len)); //base of 2 for number
	int minSize = target_len;
	int maxSize = pow(2, base + 1);
	//	clock_t tenth_of_sec_clock = CLOCKS_PER_SEC*0.01;

	vector<int> possibleValues, powers;
	double min_op;
	bool good_int;

	powers.resize(ListOfAllowedPrimes.size());
	
	//the following rules have been found to ensure that the a given length has been optimized for speed by fftw
	//1) the number must be divisible by 2
	//2) the number cannot be divisible by prime numbers other than 2,3,5and 7
	//3) 2 must be the prime number with the greatest number of factors: (x in 2^x must be greater than y in 3^y)	
	for (int i = minSize; i <= maxSize; i++){ //find the maximum possible hit
		good_int = IsDivisible(i, ListOfAllowedPrimes, powers);
		if (powers[0] >= *max_element(powers.begin(), powers.end()) && good_int){
			possibleValues.push_back(i);
			if (i / pow(2, powers[0]) < 9){ //if the remainder after dividing all powers of 2 from n is less than 9, then this has always been seen to be a very good vector length for FFT program, that is, there is no number greater than this number that results in a faster computation time				
				maxSize = i;
				break;
			}
		}
	}
	
	vector<vector<double>> numPlanFLOPS(possibleValues.size(), vector<double>(2));
	if (possibleValues.size() == 1){
		maxSize = possibleValues[0];
		min_op = 0;
	}
	else{
		for (int i = 0; i< possibleValues.size(); i++){
			numPlanFLOPS[i][0] = possibleValues[i];
			numPlanFLOPS[i][1] = AnalyzeFFTTime(possibleValues[i]);
		}

		sort(numPlanFLOPS.begin(), numPlanFLOPS.end(), [](const vector<double> & a, const vector<double> & b){ return (a[1] < b[1]); });
		maxSize = numPlanFLOPS[0][0];
		min_op = numPlanFLOPS[0][1];
	}
		
	return maxSize;

}





//model 1= sequence model
//model 2 = germline model
/*%this functino assumes => NO ALIASING BETWEEN TO SIGNALS...THAT WAY WE KNOW
%ANY OVERLAP IS ONLY DUE TO A SPECIFIC ALIGNMENT.  Any signals that overlap
%with one another would result in us not really knowing whcih overlapping sequenes
%are truly leading to the correct overlap
%i.e. compare these signal pairs [0,0,0,1,2,3];[1,2,3,0,0,0]
%WITH [1,2,3];[1,2,3];
%model1 => this defines how the nucleotide sequence that is not shifting is
%currently aligned.  it lists all the positions where a dna sequence
%exists.  For example if this sequence is put first as [ACTG.....] WHERE
%"." represents no sequence or "0 padding" then model1 would be
%[1,2,3,4,0,0,0,0,0] where 4 represents position 4 or nucleotide 4 along
%sequence
%model 2=> this defines how the nucleotide sequence that is shifting along
%the first one is aligned/designed.  For example if model 2 is [.....ACTG]
%where '." represents regions of no sequence or "0" padding, then we right
%it as [0,0,0,0,0,1,2,3,4] */
void ffthelper::MapCrossCorrelationToShift(int **&model, int sizeModel, const vector<int> & seq1, const vector<int> & seq2)
{
	vector<int> model2shift(sizeModel);

	/*if alignment Model is not empty, then go through it and remove anymemory allocated to it*/
	if (model != NULL){
		int i = 0;
		while (model[i] != NULL){
			std::free(model[i]);
			i += 1;
		}
		std::free(model);
		model = NULL;
	}

	model = (int**)malloc((sizeModel + 1)*sizeof(int*));
	for (int i = 0; i < sizeModel; i++)
		model[i] = (int*)malloc(6 * sizeof(int));

	model[sizeModel] = NULL;

	int shiftPos;
	vector<int> algn(sizeModel);

	int currentPos = 1;
	int start = -1;
	int maxLen = 0;
	int maxP = -1;
	int endP;
	int batch = -1;
	vector<int> temp(6);
	vector<vector<int> > posS;
	vector<int> temp2(3);
	int skip;
	int algnIndLeft, algnIndRight;
	for (int shiftI = 0; shiftI < sizeModel; ++shiftI){
		for (int i = 0; i < sizeModel; ++i){
			shiftPos = ((i - shiftI) + sizeModel) % sizeModel;
			model2shift[i] = seq2[shiftPos];
		}
		for (int i = 0; i < sizeModel; ++i){
			algn[i] = seq1[i] * model2shift[i];
		}
		start = -1;
		endP = -1;
		maxP = -1;
		maxLen = 0;
		skip = 0;
		for (int s = 0; s < sizeModel; ++s){
			if (algn[s] != 0){
				if (start == -1)
					start = s;
			}
			else{
				if (start != -1){
					endP = s - 1;
					batch = batch + 1;
					temp2[0] = start;
					temp2[1] = endP;
					temp2[2] = endP - start + 1;
					posS.push_back(temp2);
					if (temp2[2] > maxLen){
						maxP = batch;
						maxLen = temp2[2];
					}
					start = -1;
					endP = -1;
				}
			}
		}
		if (algn[sizeModel - 1] != 0){
			if (start == -1){
				start = sizeModel - 1;
				endP = sizeModel - 1;
			}
			else{
				endP = sizeModel - 1;
			}
			batch = batch + 1;
			temp2[0] = start;
			temp2[1] = endP;
			temp2[2] = endP - start + 1;
			posS.push_back(temp2);
			if (temp2[2] > maxLen)
				maxP = batch;
		}
		temp[0] = shiftI;
		if (maxP == -1){
			temp[1] = 0;
			temp[2] = 0;
			temp[3] = 0;
			temp[4] = 0;
			temp[5] = 0;
		}
		else{
			algnIndLeft = posS[maxP][0];
			algnIndRight = posS[maxP][1];
			temp[1] = seq1[algnIndLeft];
			temp[2] = seq1[algnIndRight];
			temp[3] = model2shift[algnIndLeft];
			temp[4] = model2shift[algnIndRight];
			temp[5] = temp[4] - temp[3] + 1;
		}
		
		for (int h = 0; h < 6; h++)
			model[shiftI][h] = temp[h];		
	}
}

//inline void ffthelper::CrossCorrelation(structvars::GermlineAlignmentSettings & alignment,  complex_t* seqSignal, complex_t* db_signal){
//void algnfunct::MakeSeqPeptide(int *list_of_peptides, int *algn_pos_indices, const int &  germ_cluster_index, complex_t *seq_complex, structvars::PeptideFFT & seq_peptide_data){
//	int algn_len;
//	//algn_pos_indices[0] = front_pos(1,1)
//	//algn_pos_indices[1] = front_pos(1,2)
//	//algn_pos_indices[2] = germPos(1,1)
//	//algn_pos_indices[3] = diagonal (germ-start)
//	complex_t *query_peptide, *germ_peptide;
//	if (list_of_peptides[algn_pos_indices[0]] == 0){ //the FFT signal of this peptide has not been defined yet
//		//algn_len = seq_peptide_data.peptide_len + 2 * (seq_peptide_data.length_flank_gap + seq_peptide_data.max_gap_length);// algn_pos_indices[1] - algn_pos_indices[0] + 1;
//		algn_len = algn_pos_indices[1] - algn_pos_indices[0] + 1;
//		seq_peptide_data.ConvertSeqFFT(seq_complex, seq_peptide_data.sequence[algn_pos_indices[0]], algn_pos_indices[0], algn_len);
//		seq_peptide_data.seq_pos_store[algn_pos_indices[0]][0] = algn_pos_indices[0];
//		//seq_peptide_data.seq_pos_store[algn_pos_indices[0]][1] = algn_pos_indices[0][1];
//		for (int i = 0; i <= seq_peptide_data.length_flank_gap; i++)
//			list_of_peptides[i + algn_pos_indices[0]] = algn_pos_indices[0] + 1;
//	}
//	algn_pos_indices[0] = list_of_peptides[algn_pos_indices[0]] - 1;
//	seq_peptide_data.queryPeptide = seq_peptide_data.sequence[algn_pos_indices[0]];
//
//	seq_peptide_data.germPeptide = seq_peptide_data.germlines[germ_cluster_index][algn_pos_indices[2]];
//	seq_peptide_data.expectedPeak = (algn_pos_indices[2] - algn_pos_indices[0]) - algn_pos_indices[3];
//}
//
//bool algnfunct::PeptideFFT(structvars::PeptideFFT & seq_peptide_data, int *diags, const double & stdevAlgn){
//	ffthelper::CrossCorrelation(seq_peptide_data.crossCorrVector, seq_peptide_data.germPeptide, seq_peptide_data.queryPeptide, seq_peptide_data.true_algn_len);
//	fftw_execute(seq_peptide_data.pep_r1);
//	double cutoff = 3 * stdevAlgn;
//	bool found = false;
//	int i = -1 * seq_peptide_data.max_gap_length;
//	int shift;
//	while (*(diags + 0) == 0 && i < 0){
//		shift = (seq_peptide_data.true_algn_len + seq_peptide_data.expectedPeak + i) % seq_peptide_data.true_algn_len;
//		if (round(seq_peptide_data.alignmentVector[shift][0])>cutoff){
//			*(diags + 0) = -1 * i; //always want diagonals to be positive
//			found = true;
//		}
//		i++;
//	}
//
//	if (!found){
//		i = seq_peptide_data.max_gap_length;
//		while (*(diags + 1) == 0 && i > 0){
//			shift = (seq_peptide_data.true_algn_len + seq_peptide_data.expectedPeak + i) % seq_peptide_data.true_algn_len;
//			if (round(seq_peptide_data.alignmentVector[shift][0]) > cutoff){
//				*(diags + 1) = i;
//				found = true;
//			}
//			i--;
//		}
//	}
//
//	if (!found)
//		found = round(seq_peptide_data.alignmentVector[seq_peptide_data.expectedPeak][0]) > cutoff ? true : false;
//	return found;
//}