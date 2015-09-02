#include "FFTAlign.h"

using namespace std;

FFTAlign::FFTAlign(int n){
	alignmentLength = n; //this is the desired length for performing fft (usually its defined as max seq len + max germ len)
	/****initalize the complex vars/pointers using NULL (this is how we will now if they are being created for teh first time to avoid segmentation error/memory leaks)***/
	forwardT = NULL;
	reverseT = NULL;
	inputVector = NULL;
	fourierVector = NULL;
	crossCorrVector = NULL;
	alignmentVector = NULL;
	expectedErrorModel = NULL;
	alignmentModel =  NULL;
	/*****************/

	vector<int> seqModel(n), germModel(n);
	for (int i = 0; i < n; i++){
		seqModel[i] = i;
		germModel[i] = i;
	}

	UpdateOverlapModel(seqModel, germModel);
	CreateFFTPlan();
	ModelBinomialAlignmentError();
}

FFTAlign::FFTAlign(int n, structvars::FFTSettings alignment_settings){
	alignmentLength = n; //this is the desired length for performing fft (usually its defined as max seq len + max germ len)
	fftparams = alignment_settings;

	/****initalize the complex vars/pointers using NULL (this is how we will now if they are being created for teh first time to avoid segmentation error/memory leaks)***/
	forwardT = NULL;
	reverseT = NULL;
	inputVector = NULL;
	fourierVector = NULL;
	crossCorrVector = NULL;
	alignmentVector = NULL;
	expectedErrorModel = NULL;
	alignmentModel = NULL;
	/*****************/
	vector<int> seqModel(n), germModel(n);
	
	for (int i = 0; i < n; i++){
		seqModel[i] = i;
		germModel[i] = i;
	}

	UpdateOverlapModel(seqModel, germModel);
	CreateFFTPlan();
	ModelBinomialAlignmentError();
}

FFTAlign::FFTAlign(int n, const std::vector<int> & seqModel, const std::vector<int> & germModel, structvars::FFTSettings alignment_settings){
	alignmentLength = n; //this is the desired length for performing fft (usually its defined as max seq len + max germ len)
	fftparams = alignment_settings;

	/****initalize the complex vars/pointers using NULL (this is how we will now if they are being created for teh first time to avoid segmentation error/memory leaks)***/
	forwardT = NULL;
	reverseT = NULL;
	inputVector = NULL;
	fourierVector = NULL;
	crossCorrVector = NULL;
	alignmentVector = NULL;
	alignmentModel = NULL;
	expectedErrorModel = NULL;
	/*****************/

	UpdateOverlapModel(seqModel, germModel);
	CreateFFTPlan();
	ModelBinomialAlignmentError();
}

FFTAlign::FFTAlign(int n, const std::vector<int> & seqModel, const std::vector<int> & germModel){
	alignmentLength = n; //this is the desired length for performing fft (usually its defined as max seq len + max germ len)	

	/****initalize the complex vars/pointers using NULL (this is how we will now if they are being created for teh first time to avoid segmentation error/memory leaks)***/
	forwardT = NULL;
	reverseT = NULL;
	inputVector = NULL;
	fourierVector = NULL;
	crossCorrVector = NULL;
	alignmentVector = NULL;
	alignmentModel = NULL;
	expectedErrorModel = NULL;
	/*****************/

	UpdateOverlapModel(seqModel, germModel);
	CreateFFTPlan();
	ModelBinomialAlignmentError();
}

FFTAlign::FFTAlign(const FFTAlign & fftcopy){

	alignmentLength = fftcopy.alignmentLength;

	fftparams = fftcopy.fftparams;
	crossCorrVector = (complex_t*)malloc(alignmentLength*sizeof(complex_t));
	inputVector = (complex_t*)malloc(alignmentLength*sizeof(complex_t));
	fourierVector = (complex_t*)malloc(alignmentLength*sizeof(complex_t));
	alignmentVector = (complex_t*)malloc(alignmentLength*sizeof(complex_t));

	memcpy(crossCorrVector, fftcopy.crossCorrVector, alignmentLength*sizeof(complex_t));
	memcpy(inputVector, fftcopy.inputVector, alignmentLength*sizeof(complex_t));
	memcpy(fourierVector, fftcopy.fourierVector, alignmentLength*sizeof(complex_t));
	memcpy(alignmentVector, fftcopy.alignmentVector, alignmentLength*sizeof(complex_t));

	//make fft plans
	forwardT = fftw_plan_dft_1d(alignmentLength, inputVector, fourierVector, FFTW_FORWARD, FFTW_PATIENT); //forward transform
	reverseT = fftw_plan_dft_1d(alignmentLength, crossCorrVector, alignmentVector, FFTW_BACKWARD, FFTW_PATIENT); //reverse fft TRANSFORM

	expectedErrorModel = (double**)malloc(3 * sizeof(double*));
	expectedErrorModel = fftcopy.expectedErrorModel;
	//memcpy(expectedErrorModel, fftcopy.expectedErrorModel, 3 * sizeof(double*));
	//for (int i = 0; i < 3; i++){
	//	expectedErrorModel[i] = (double*)malloc((alignmentLength + 1) * sizeof(double));
	//	memcpy(expectedErrorModel[i], fftcopy.expectedErrorModel[i], (alignmentLength + 1) * sizeof(double));
	//}

	alignmentModel = (int**)malloc(alignmentLength*sizeof(int*));
	alignmentModel = fftcopy.alignmentModel;
	
	//for (int i = 0; i < alignmentLength; i++)
		//memcpy(alignmentModel[i], fftcopy.alignmentModel[i], 6 * sizeof(int));
}

FFTAlign::~FFTAlign(){
	if (inputVector != NULL){
		std::free(inputVector);
		inputVector = NULL;
	}
	if (fourierVector != NULL){
		std::free(fourierVector);
		fourierVector = NULL;
	}
	if (crossCorrVector != NULL){
		std::free(crossCorrVector);
		crossCorrVector = NULL;
	}
	if (alignmentVector != NULL){
		std::free(alignmentVector);
		alignmentVector = NULL;
	}

	if (alignmentModel != NULL){
		for (int i = 0; i < alignmentLength; i++){
			std::free(alignmentModel[i]);
		}				
		std::free(alignmentModel);
	}
	
	fftw_destroy_plan(forwardT);
	fftw_destroy_plan(reverseT);

	for (int i = 0; i < 3; i++)
		std::free(expectedErrorModel[i]);
	std::free(expectedErrorModel);
	expectedErrorModel = NULL;
}

/*PUBLIC FUNCTIONS*/
int FFTAlign::GetAlignmentLength(){
	return alignmentLength;
}

//same fucntion as above but allowes user to send in a new model manually
void FFTAlign::UpdateOverlapModel(const vector<int> & seqModel, const vector<int> & germModel){
	vector<int> seqModel2 = seqModel, germModel2 = germModel;
	
	assert(seqModel2.size() == germModel2.size());
	if (seqModel2.size() != alignmentLength){
		printf("\nWARNING: THE MODEL PROVIDED TO THE FFT ALIGN CLASS IS NOT THE CORRECT LENGTH, RESIZING...THIS MAY CAUSE PROBLEMS\n");
		seqModel2.resize(alignmentLength);
		germModel2.resize(alignmentLength);
	}
		
	ffthelper::MapCrossCorrelationToShift(alignmentModel, alignmentLength, seqModel2, germModel2);
}

int FFTAlign::AlignSequences(complex_t *sequence, complex_t *germline){ //actually performs an fft alignment
	ffthelper::CrossCorrelation(crossCorrVector, germline, sequence, alignmentLength);
	fftw_execute(reverseT);
	
	bestDiag = 0;
	bestDiagScore = 0;

	for (int i = 0; i < alignmentLength; i++){
		alignmentVector[i][0] = std::round(alignmentVector[i][0]);
		if (alignmentVector[i][0]>bestDiagScore){
			bestDiagScore = alignmentVector[i][0];
			bestDiag = i;
		}
	}

	//erase any of these variables from previously alignemnts/functionc alls
	totalScore = bestDiagScore;
	gappedScore = bestDiagScore;
	gapPenalty = 0;
	minDiag = 0;
	maxDiag = 0;
	numGaps = 0;
	return bestDiag;
}

//look for other non random alignemnts/diagonals with strong alignments to the sequence.
//presence of mulitple alignments indicates gaps in the sequence
//alignmentModel maps the position/index of a cross correlation to the alignment length between sequences (col_num => the column number which contains the alignment length info)
void FFTAlign::IdentifyAdditionalAlignments(int max_len, int col_num){
	minDiag = 0;
	maxDiag = 0;
	gapPenalty = 0;
	totalScore = 0;
	numGaps = 0;
	gappedScore = 0;

	int newgap = 0;
	int prevGap = -1;
	int algn_len;
	int minAlgnPosition = max(bestDiag - fftparams.maxGap, 0); //look at scroes MAXGAP position below the best diagonal identified
	int maxAlgnPosition = min(bestDiag + fftparams.maxGap, alignmentLength - 1); ////look at scroes MAXGAP position above the best diagonal identified
	double threshold;
	for (int k = minAlgnPosition; k <= maxAlgnPosition; k++){
		algn_len = min(alignmentModel[k][col_num], max_len);
		threshold = expectedErrorModel[2][algn_len]; //alignmentModel[k][5] => position/index of k and column 5 represents the alignment length at that position
		if (alignmentVector[k][0] > threshold){ //then we found a hit
			totalScore += alignmentVector[k][0];
			if (prevGap == -1)
				minDiag = k;
			else{
				newgap = k - prevGap;
				numGaps += newgap;
				gapPenalty += fftparams.fftGapOpen + fftparams.fftExtendGap*(newgap - 1);
			}
			maxDiag = k;
			prevGap = k;
		}
	}

	minDiag = bestDiag - minDiag;
	maxDiag = maxDiag - bestDiag;
	gappedScore = totalScore - gapPenalty;
}

/*ACCESSORS*/
double FFTAlign::ReturnBinomialError(int errorType, int length){ //first int is the row number, second int is the column number (alignment length in this case)
	return expectedErrorModel[errorType][length];
}

double FFTAlign::ReturnAlignmentScoreAtPosition(int diagPosition){
	return alignmentVector[diagPosition][0]; //only care about real numbers (not imaginary)
}

/*PRIVATE FUNCTIONS*/
//Initializing/making an FFTW plan for the alignment.
void FFTAlign::CreateFFTPlan(){
	if (forwardT != NULL) //this is a safety check, just incase we have already made a plan
		fftw_destroy_plan(forwardT);
	if (reverseT != NULL)
		fftw_destroy_plan(reverseT);
	if (inputVector != NULL)
		std::free(inputVector);
	if (fourierVector != NULL)
		std::free(fourierVector);

	if (crossCorrVector != NULL)
		std::free(crossCorrVector);
	if (alignmentVector != NULL)
		std::free(alignmentVector);
	
	/*if (!findIdealAlgnSize)
		alignmentLength = minAlignmentLength;
	else{
		alignmentLength = ffthelper::FindIdealAlgnSize(minAlignmentLength); //the given alignment length is probably not ideal speed for FFT, so we will look for a nearby length
		//alignmentLength = minAlignmentLength;
		//printf("DONT FORGET TO CHANGE ALIGNMENT LENGTH SETTING IN CREATE FFT PLAN");
	}*/
	

	//allocate the proper memory before making fft plans
	inputVector = (complex_t*)malloc(alignmentLength*sizeof(complex_t));
	fourierVector = (complex_t*)malloc(alignmentLength*sizeof(complex_t));
	crossCorrVector = (complex_t*)malloc(alignmentLength*sizeof(complex_t));
	alignmentVector = (complex_t*)malloc(alignmentLength*sizeof(complex_t));

	//make fft plans
	forwardT = fftw_plan_dft_1d(alignmentLength, inputVector, fourierVector, FFTW_FORWARD, FFTW_PATIENT); //forward transform
	reverseT = fftw_plan_dft_1d(alignmentLength, crossCorrVector, alignmentVector, FFTW_BACKWARD, FFTW_PATIENT); //reverse fft TRANSFORM
	return;
}

//We ue this function to estimate the background error of the FFT alignment
//For now we do not allow the user to change the "mathc and mismatch" parameters for FFT alignments
void FFTAlign::ModelBinomialAlignmentError(){
	double FFTMatchScore = 1.000, FFTMismatchScore = -1.000;
	double varianceMatch = pow(FFTMatchScore, 2)*0.25*(1 - 0.25);
	double varianceMisMatch = pow(FFTMismatchScore, 2)*0.25*(1 - 0.25);
	double covariance = -2 * (0.25)*(0.25);
	double expScoreDeviation = sqrt(1 * (varianceMatch + varianceMisMatch - covariance));

	/*remove any preallocated memory for a previous error model*/
	if (expectedErrorModel != NULL){
		std::free(expectedErrorModel[0]); //row 0 is the expected score
		std::free(expectedErrorModel[1]); //row 1 is the expected score + variance
		std::free(expectedErrorModel[2]); //row 2 is column 1*the sensitivity
	}

	/*make an error model. this will predict the expected variance of rnadom alignmetn between DNA sequences of specified length*/
	expectedErrorModel = (double **)malloc((3)*sizeof(double*));
	expectedErrorModel[0] = (double*)malloc((alignmentLength + 1)*sizeof(double));
	expectedErrorModel[1] = (double*)malloc((alignmentLength + 1)*sizeof(double));
	expectedErrorModel[2] = (double*)malloc((alignmentLength + 1)*sizeof(double));

	for (int i = 0; i <= alignmentLength; i++){
		expectedErrorModel[0][i] = (1.0000*expScoreDeviation*sqrt(i)); //at i = 0, length is 1, not 0 //standard error in alignment
		expectedErrorModel[1][i] = (0.2500*(i)+sqrt(i)*sqrt(varianceMatch)); //mean + standard error in alignment
		expectedErrorModel[2][i] = fftparams.sensitivity*expectedErrorModel[0][i];//cutoff for identifying non-random alignments
	}

	return;
}

bool FFTAlign::FindGappedPeaks(int * info, double sensitivity){
	int algn_len = info[2];
	int expectedPeak = info[3];
	double cutoff =  sensitivity * expectedErrorModel[0][algn_len];
	bool found = false;

	int i = -1 * fftparams.maxGap;
	int shift;

	info[0] = 0;//mindiag
	info[1] = 0;//maxdiag

	while (info[0] == 0 && i < 0){
		shift = (alignmentLength + expectedPeak + i) % (alignmentLength);

		if (alignmentVector[shift][0]>cutoff){
			info[0] = -1 * i; //always want diagonals to be positive
			found = true;
		}
		i++;
	}

	//if (!found){
		i = fftparams.maxGap;
		while (info[1] == 0 && i > 0){
			shift = (alignmentLength + expectedPeak + i) % alignmentLength;
			if (alignmentVector[shift][0] > cutoff){
				info[1] = i;
				found = true;
			}
			i--;
		}
	//}

	if (!found)
		found = round(alignmentVector[expectedPeak][0]) > cutoff ? true : false;
	return found;
}