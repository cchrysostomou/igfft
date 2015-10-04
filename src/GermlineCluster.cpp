#include "GermlineCluster.h"

using namespace::std;

GermlineCluster::GermlineCluster(){
	clusterIndex = -1;
	consensus = "";
	paddedNTSequence = "";
	numGermlineMembers = 0;
	clusterSeqLength = 0;
	querySeqToBeModified = "";
	minPerIdenPetideAlign = 0.8;

	complexSeqPad = NULL;
	fourierSeq = NULL;
	germlineIndices = NULL;
	germlineStart = NULL;
	germlineSequences = NULL;
	peptideFFT = NULL;
	sw_align_cluster = NULL;
	queryInfoPointer = NULL;
	seqIns = NULL;
	seqDel = NULL;

	numPeptides = 0;
	minDiag = 0;
	maxDiag = 0;
	numFFTGaps = 0;
	FFTgappedScore = 0;
	FFTgaplessScore = 0;
	SWScore = 0;
	finalScore = 0;
	numSteps = 0;
	lineExit = 0;
	ds_us_distance = 50;
	numFFTPlans = 0;
	for (int i = 0; i < 5; i++)
		alignToSegment[i] = true;
}

GermlineCluster::GermlineCluster(int cluster, const std::string & c, const structvars::AlignmentProgramSettings & params){
	clusterIndex = cluster;
	consensus = c;
	numGermlineMembers = 0;
	clusterSeqLength = c.length();
	paddedNTSequence = "";
	querySeqToBeModified = "";
	algnParams = params;
	minPerIdenPetideAlign = 0.8;

	complexSeqPad = NULL;
	fourierSeq = NULL;
	germlineIndices = NULL;
	germlineStart = NULL;
	germlineSequences = NULL;
	peptideFFT = NULL;
	sw_align_cluster = NULL;
	queryInfoPointer = NULL;
	seqIns = NULL;
	seqDel = NULL;

	numPeptides = 0;
	minDiag = 0;
	maxDiag = 0;
	numFFTGaps = 0;
	FFTgappedScore = 0;
	FFTgaplessScore = 0;
	SWScore = 0;
	finalScore = 0;
	numSteps = 0;
	lineExit = 0;
	ds_us_distance = 50;
	numFFTPlans = 0;
	for (int i = 0; i < 5; i++)
		alignToSegment[i] = true;
}

GermlineCluster::GermlineCluster(int cluster, const std::string & c, const structvars::AlignmentProgramSettings & params, const std::vector<int> & germlineMembers, const std::vector<structvars::Germline> & germDB, const std::vector<structvars::PairwiseClusterAlignment> & alignment_to_clusters){
	clusterIndex = cluster;
	consensus = c;
	numGermlineMembers = germlineMembers.size();
	alignmentToOtherClusters = alignment_to_clusters;
	clusterSeqLength = c.length();
	paddedNTSequence = "";
	querySeqToBeModified = "";
	algnParams = params;
	minPerIdenPetideAlign = 0.8;

	germlineSequences = new string[numGermlineMembers];
	germlineStart = new int[numGermlineMembers];
	germlineIndices = new int[numGermlineMembers];

	for (int i = 0; i < numGermlineMembers; i++){
		germlineIndices[i] = germlineMembers[i];
		germlineSequences[i] = germDB[germlineIndices[i]].germlineSeqData.seq;
		germlineStart[i] = germDB[germlineIndices[i]].seqStart;
	}

	complexSeqPad = NULL;
	fourierSeq = NULL;
	peptideFFT = NULL;
	sw_align_cluster = NULL;
	queryInfoPointer = NULL;
	seqIns = NULL;
	seqDel = NULL;

	numPeptides = 0;
	minDiag = 0;
	maxDiag = 0;
	numFFTGaps = 0;
	FFTgappedScore = 0;
	FFTgaplessScore = 0;
	SWScore = 0;
	finalScore = 0;
	numSteps = 0;
	lineExit = 0;
	ds_us_distance = 50;
	numFFTPlans = 0;
	for (int i = 0; i < 5; i++)
		alignToSegment[i] = true;
}

GermlineCluster::~GermlineCluster(){
	if (germlineIndices != NULL)
		delete[] germlineIndices;
	if (germlineStart != NULL)
		delete[] germlineStart;
	if (germlineSequences != NULL)
		delete[] germlineSequences;
	if (complexSeqPad != NULL)
		std::free(complexSeqPad);
	if (fourierSeq != NULL){
		for (int i = 0; i < numFFTPlans;i++)
			std::free(fourierSeq[i]);
		std::free(fourierSeq);
	}
	if (peptideFFT != NULL){
		assert(numPeptides > 0);
		for (int i = 0; i < numPeptides; i++)
			std::free(peptideFFT[i]);
		std::free(peptideFFT);
	}

	if (seqIns != NULL){
		for (int k = 0; k < 5 * algnParams.swParams.maxGap; k++)
			delete[] seqIns[k];
		delete[] seqIns;
	}

	if (seqDel != NULL){
		for (int k = 0; k < 5 * algnParams.swParams.maxGap; k++)
			delete[] seqDel[k];
		delete[] seqDel;
	}

	if (sw_align_cluster != NULL)
		delete(sw_align_cluster);
}

void GermlineCluster::InitializeComplex(const vector<int> & algnLens){
	if (complexSeqPad != NULL)
		std::free(complexSeqPad);

	if (fourierSeq != NULL){
		for (int i = 0; i < numFFTPlans; i++)
			std::free(fourierSeq[i]);
		std::free(fourierSeq);
	}
	numFFTPlans = algnLens.size();

	complexSeqPad = (complex_t*)malloc(algnLens[numFFTPlans-1]*sizeof(complex_t));
	
	fourierSeq = (complex_t**)malloc(numFFTPlans*sizeof(complex_t*));
	for (int i = 0; i<numFFTPlans;i++)
		fourierSeq[i] = (complex_t*)malloc(algnLens[numFFTPlans - 1] * sizeof(complex_t));
}

/* PRIVATE MODIFIERS*/
void GermlineCluster::AddSWAlignmentInfo(int maxQueryLength, int maxClusterLength){
	sw_align_cluster = new SWAlignment(maxQueryLength, maxClusterLength, algnParams.swParams);
}

void GermlineCluster::AddGermlineMemberInfo(const std::vector<int> & germlineMembers, const std::vector<structvars::Germline> & germDB){
	if (germlineIndices != NULL)
		delete(germlineIndices);
	if (germlineStart != NULL)
		delete(germlineStart);
	if (germlineSequences != NULL)
		delete(germlineSequences);

	numGermlineMembers = germlineMembers.size();
	germlineSequences = new string[numGermlineMembers];
	germlineStart = new int[numGermlineMembers];
	germlineIndices = new int[numGermlineMembers];
	locusName = germDB[germlineMembers[0]].locus;
	for (int i = 0; i < numGermlineMembers; i++){
		germlineIndices[i] = germlineMembers[i];
		germlineSequences[i] = germDB[germlineIndices[i]].germlineSeqData.seq;
		germlineStart[i] = germDB[germlineIndices[i]].seqStart;
		if (germDB[germlineIndices[i]].locus != locusName){
			printf("Multiple loci have been grouped into the same germline clusters\n");
		}
	}
}

void GermlineCluster::MakePeptideFFT(int maxGermLen){
	if (peptideFFT != NULL){
		assert(numPeptides > 0);
		for (int i = 0; i < numPeptides; i++)
			std::free(peptideFFT[i]);
		std::free(peptideFFT);
	}

	complex_t *tempPeptide = (complex_t*)malloc((*pointer_to_peptide_alignment).GetAlignmentLength()*sizeof(complex_t));
	numPeptides = maxGermLen;// -algnParams.peptideParams.peptide_len + 1;
	peptideFFT = (complex_t**)malloc(numPeptides*sizeof(complex_t*));
	int pepLen;
	for (int i = 0; i < numPeptides; i++){
		pepLen = i + algnParams.peptideParams.peptide_len > maxGermLen ? maxGermLen - i : algnParams.peptideParams.peptide_len;
		memcpy(tempPeptide, complexSeqPad + i, pepLen*sizeof(complex_t));
		for (int k = pepLen; k < (*pointer_to_peptide_alignment).GetAlignmentLength(); k++){
			tempPeptide[k][0] = 0;
			tempPeptide[k][1] = 0;
		}
		peptideFFT[i] = (complex_t*)malloc((*pointer_to_peptide_alignment).GetAlignmentLength()*sizeof(complex_t));
		(*pointer_to_peptide_alignment).ForwardFFT(tempPeptide, peptideFFT[i]);
		for (int k = 0; k < (*pointer_to_peptide_alignment).GetAlignmentLength(); k++){
			peptideFFT[i][k][0] /= (*pointer_to_peptide_alignment).GetAlignmentLength();
			peptideFFT[i][k][1] /= (-1 * (*pointer_to_peptide_alignment).GetAlignmentLength());
		}
	}
	std::free(tempPeptide);
}

void GermlineCluster::AddPairwiseClusterInfo(const std::vector<structvars::PairwiseClusterAlignment> & alignment_to_clusters){
	alignmentToOtherClusters = alignment_to_clusters;
}

bool GermlineCluster::QuickAlignCheck(region r){
	//if (!alignToSegment[r])
	//	return true;

	int algnInfo[3];
	int scores[2] = { 1, 0 };
	switch (r){
	case BEG:
		algnInfo[0] = alignmentToSequence[1];
		algnInfo[1] = alignmentToSequence[3];
		algnInfo[2] = algnParams.peptideParams.peptide_len;
		break;
	case BEG_DS:
		algnInfo[0] = alignmentToSequence[1] + ds_us_distance - 1;
		algnInfo[1] = alignmentToSequence[3] + ds_us_distance - 1;
		algnInfo[2] = algnParams.peptideParams.peptide_len;
		break;
	case MID:
		algnInfo[0] = (alignmentToSequence[1] + alignmentToSequence[5] / 2);
		algnInfo[1] = (alignmentToSequence[3] + alignmentToSequence[5] / 2);
		algnInfo[2] = algnParams.peptideParams.peptide_len;
		break;
	case US_STOP:
		algnInfo[0] = alignmentToSequence[2] - algnParams.peptideParams.peptide_len - ds_us_distance + 1;
		algnInfo[1] = alignmentToSequence[4] - algnParams.peptideParams.peptide_len - ds_us_distance + 1;
		algnInfo[2] = algnParams.peptideParams.peptide_len;

		break;
	case LAST:
		algnInfo[0] = alignmentToSequence[2] - algnParams.peptideParams.peptide_len + 1;
		algnInfo[1] = alignmentToSequence[4] - algnParams.peptideParams.peptide_len + 1;
		algnInfo[2] = algnParams.peptideParams.peptide_len;
		break;
	default:
		return true;
	}
	double algnScore = dnafunctions::AlignDiagonal(algnInfo, querySeqToBeModified, paddedNTSequence, scores);
	//length of alignment is changed in AlignDiagonal function. this is because seomtimes alignment length isnot always "peptide_len"
	return algnScore >= minPerIdenPetideAlign*algnInfo[2];// (*pointer_to_full_seq_alignment).expectedErrorModel[1][algnInfo[2]];
}

bool GermlineCluster::QuickAlignCheck(region r, const GermlineCluster & refSeq){
	int algnInfo[3];
	int scores[2] = { 1, 0 };
	switch (r){
	case BEG:
		algnInfo[0] = refSeq.alignmentToSequence[1];
		algnInfo[1] = refSeq.alignmentToSequence[3];
		algnInfo[2] = refSeq.algnParams.peptideParams.peptide_len;
		break;
	case BEG_DS:
		algnInfo[0] = refSeq.alignmentToSequence[1] + ds_us_distance - 1;
		algnInfo[1] = refSeq.alignmentToSequence[1] + ds_us_distance - 1;
		algnInfo[2] = refSeq.algnParams.peptideParams.peptide_len;
		break;
	case MID:
		algnInfo[0] = (refSeq.alignmentToSequence[1] + refSeq.alignmentToSequence[5] / 2);
		algnInfo[1] = (refSeq.alignmentToSequence[3] + refSeq.alignmentToSequence[5] / 2);
		algnInfo[2] = algnParams.peptideParams.peptide_len;
		break;
	case US_STOP:
		algnInfo[0] = refSeq.alignmentToSequence[2] - refSeq.algnParams.peptideParams.peptide_len - ds_us_distance + 1;
		algnInfo[1] = refSeq.alignmentToSequence[4] - refSeq.algnParams.peptideParams.peptide_len - ds_us_distance + 1;
		algnInfo[2] = algnParams.peptideParams.peptide_len;

		break;
	case LAST:
		algnInfo[0] = refSeq.alignmentToSequence[2] - refSeq.algnParams.peptideParams.peptide_len + 1;
		algnInfo[1] = refSeq.alignmentToSequence[4] - refSeq.algnParams.peptideParams.peptide_len + 1;
		algnInfo[2] = algnParams.peptideParams.peptide_len;
		break;
	default:
		return true;
	}
	double algnScore = dnafunctions::AlignDiagonal(algnInfo, querySeqToBeModified, paddedNTSequence, scores);
	//length of alignment is changed in AlignDiagonal function. this is because seomtimes alignment length isnot always "peptide_len"
	return algnScore >= minPerIdenPetideAlign*algnInfo[2];// (*pointer_to_full_seq_alignment).expectedErrorModel[1][algnInfo[2]];
}
//, FFTAlign *peptide_seq_alignment, complex_t **seqPeptideList, int *seqPosStore;
bool GermlineCluster::FindDiagInSpecificRegion(region r){
	if (!alignToSegment[r])
		return true;
	int seqPos = algnPepResults[r][QUERY][FIRST];
	int germPos = algnPepResults[r][GERMLINE][FIRST];
	int algn_len = algnPepResults[r][QUERY][FINAL] - algnPepResults[r][QUERY][FIRST] + 1;
	if ((*pointer_to_seq_peptide_list)[seqPos] == 0){
		complex_t *temp = (complex_t*)malloc((*pointer_to_peptide_alignment).alignmentLength*sizeof(complex_t));
		memcpy(temp, (*queryInfoPointer).complSeq + seqPos, algn_len*sizeof(complex_t));
		for (int k = algn_len; k < (*pointer_to_peptide_alignment).alignmentLength; k++){
			temp[k][0] = 0;
			temp[k][1] = 0;
		}

		(*pointer_to_peptide_alignment).ForwardFFT(temp, (*pointer_to_fft_seq_peptides)[seqPos]);
		for (int k = 0; k < algnParams.peptideParams.length_flank_gap; k++)
			(*pointer_to_seq_peptide_list)[k + seqPos] = seqPos + 1;
		std::free(temp);
	}		
	int expectedPeak = (algnPepResults[r][GERMLINE][FIRST] - ((*pointer_to_seq_peptide_list)[seqPos] - 1)) - (alignmentToSequence[3] - alignmentToSequence[1]);
	(*pointer_to_peptide_alignment).AlignSequences((*pointer_to_fft_seq_peptides)[(*pointer_to_seq_peptide_list)[seqPos] - 1], peptideFFT[germPos]);
	int diagInfo[4];
	diagInfo[0] = 0;
	diagInfo[1] = 0;
	diagInfo[2] = algnParams.peptideParams.peptide_len;//algn_len;
	diagInfo[3] = expectedPeak;
	bool foundPeak = (*pointer_to_peptide_alignment).FindGappedPeaks(diagInfo,algnParams.peptideParams.peptide_gap_sensitivity); //after aligning the peptide FFT, search the resulting alignment vector for any peaks that we assume are not due to random nucleotide alignment
	minDiag = max(minDiag, diagInfo[0]);
	maxDiag = max(maxDiag, diagInfo[1]);
	return foundPeak;
}

void GermlineCluster::DivideSeqIntoPeptides(){
	int totalGaps = algnParams.swParams.maxGap + algnParams.peptideParams.length_flank_gap;
	int min_len = 2 * totalGaps + algnParams.peptideParams.peptide_len;
	int algnLen = alignmentToSequence[2] - alignmentToSequence[1];

	//PEPTIDE ONE, WITH RESPECT TO QUERY
	int startP = alignmentToSequence[1] - totalGaps; //algnPos[1] - gapSeq
	algnPepResults[BEG][QUERY][FIRST] = startP < querySeqStart ? querySeqStart : startP;
	int stopP = algnPepResults[BEG][QUERY][FIRST] + min_len - 1; //begP(1,2) = begP(1,1)+peptideSize+2*gapSeq-1
	algnPepResults[BEG][QUERY][FINAL] = stopP>querySeqEnd - 1 ? querySeqEnd - 1 : stopP;
	//PEPTIDE ONE, WITH RESPECT TO GERMLINE
	algnPepResults[BEG][GERMLINE][FIRST] = alignmentToSequence[3];
	stopP = algnPepResults[BEG][GERMLINE][FIRST] + algnParams.peptideParams.peptide_len - 1;
	algnPepResults[BEG][GERMLINE][FINAL] = stopP >= clusterSeqLength ? clusterSeqLength - 1 : stopP;
	alignToSegment[BEG] = true; //always true

	//LAST PEPTIDE, WITH RESPECT TO QUERY
	stopP = alignmentToSequence[2] + totalGaps;
	algnPepResults[LAST][QUERY][FINAL] = stopP > querySeqEnd ? querySeqEnd : stopP;
	algnPepResults[LAST][QUERY][FIRST] = algnPepResults[LAST][QUERY][FINAL] - min_len + 1;
	//LAST PEPTIDE, WITH RESPECT TO GERMLINE
	algnPepResults[LAST][GERMLINE][FINAL] = alignmentToSequence[4];// > clusterSeqLength ? clusterSeqLength : alignmentToSequence[4];
	algnPepResults[LAST][GERMLINE][FIRST] = algnPepResults[LAST][GERMLINE][FINAL] - algnParams.peptideParams.peptide_len + 1;

	if (algnPepResults[LAST][QUERY][FIRST] <= algnPepResults[BEG][QUERY][FIRST] || algnPepResults[LAST][GERMLINE][FIRST] <= algnPepResults[BEG][GERMLINE][FIRST]){
		//if the sequence is small then no point in aliging to any more regions
		alignToSegment[LAST] = false;
		alignToSegment[BEG_DS] = false;
		alignToSegment[MID] = false;
		alignToSegment[US_STOP] = false;
		return;
	}
	else
		alignToSegment[LAST] = true;

	//MID PEPTIDE => WITH RESPECT TO QUERY
	double middleNTPos = round((algnLen / 2.0) - (algnParams.peptideParams.peptide_len / 2.0));  //middleNTPos = round(algnLen/2 - peptidesize/2)
	algnPepResults[MID][QUERY][FIRST] = algnPepResults[BEG][QUERY][FIRST] + middleNTPos - 1;
	algnPepResults[MID][QUERY][FINAL] = algnPepResults[BEG][QUERY][FINAL] + middleNTPos - 1;
	//MID PEPTIDE => WITH RESPECT TO GERMLINE
	algnPepResults[MID][GERMLINE][FIRST] = algnPepResults[BEG][GERMLINE][FIRST] + middleNTPos - 1;
	algnPepResults[MID][GERMLINE][FINAL] = algnPepResults[BEG][GERMLINE][FINAL] + middleNTPos - 1;
	if (algnPepResults[MID][QUERY][FIRST] >= algnPepResults[LAST][QUERY][FIRST] || algnPepResults[MID][GERMLINE][FIRST] >= algnPepResults[LAST][GERMLINE][FIRST]){
		//if the sequence is small then no point in aliging to any more regions
		alignToSegment[MID] = false;
		alignToSegment[US_STOP] = false;
		alignToSegment[BEG_DS] = false;
		return;
	}
	else
		alignToSegment[MID] = true;

	//DS_50 PEPTIDE => PEPTIDE TWO, WITH RESPECT TO QUERY
	algnPepResults[BEG_DS][QUERY][FIRST] = algnPepResults[BEG][QUERY][FIRST] + ds_us_distance - 1;
	algnPepResults[BEG_DS][QUERY][FINAL] = algnPepResults[BEG][QUERY][FINAL] + ds_us_distance - 1;
	//DS_50 PEPTIDE => PEPTIDE TWO, WITH RESPECT TO GERMLINE
	algnPepResults[BEG_DS][GERMLINE][FIRST] = algnPepResults[BEG][GERMLINE][FIRST] + ds_us_distance - 1;
	algnPepResults[BEG_DS][GERMLINE][FINAL] = algnPepResults[BEG][GERMLINE][FINAL] + ds_us_distance - 1;

	if (algnPepResults[BEG_DS][QUERY][FIRST] >= algnPepResults[MID][QUERY][FIRST] || algnPepResults[BEG_DS][GERMLINE][FIRST] >= algnPepResults[MID][GERMLINE][FIRST]){
		alignToSegment[BEG_DS] = false;
		alignToSegment[US_STOP] = false;
		return;
	}
	else
		alignToSegment[BEG_DS] = true;

	//US_TO PEPTIDE => PEPTIDE FOUR, WITH RESPECT TO QUERY
	algnPepResults[US_STOP][QUERY][FIRST] = algnPepResults[LAST][QUERY][FIRST] - ds_us_distance + 1;
	algnPepResults[US_STOP][QUERY][FINAL] = algnPepResults[LAST][QUERY][FINAL] - ds_us_distance + 1;
	//US_TO PEPTIDE => PEPTIDE FOUR, WITH RESPECT TO GERMLINE
	algnPepResults[US_STOP][GERMLINE][FIRST] = algnPepResults[LAST][GERMLINE][FIRST] - ds_us_distance + 1;
	algnPepResults[US_STOP][GERMLINE][FINAL] = algnPepResults[LAST][GERMLINE][FINAL] - ds_us_distance + 1;

	if (algnPepResults[US_STOP][GERMLINE][FIRST] <= algnPepResults[MID][QUERY][FIRST] || algnPepResults[US_STOP][GERMLINE][FIRST] <= algnPepResults[MID][GERMLINE][FIRST])
		alignToSegment[US_STOP] = false;
	else
		alignToSegment[US_STOP] = true;
	return;
}

void GermlineCluster::FindAllGapsInQueryToClusterAlgn(){
	numSeqIns[0] = 0;
	numSeqIns[1] = 0;
	
	numSeqDel[0] = 0;	
	numSeqDel[1] = 0;
	
	bool repeatFwd = false, repeatEnd = false, good_alignment, foundHit, problemChild=false;
	querySeqToBeModified = (*queryInfoPointer).seq;

	DivideSeqIntoPeptides();
	good_alignment = QuickAlignCheck(BEG); //v1
	//good_alignment = false; //v2
	if (good_alignment){
		repeatFwd = false;
	}
	else{
		foundHit = FindDiagInSpecificRegion(BEG);
		repeatFwd = (foundHit || !alignToSegment[BEG_DS]) ? false : true;
	}
	
	if (alignToSegment[LAST]){
		good_alignment = QuickAlignCheck(LAST);//v1
		//good_alignment = false; //v2
		if (good_alignment){
			repeatEnd = false;
		}
		else{
			foundHit = FindDiagInSpecificRegion(LAST);
			repeatEnd = (foundHit || !alignToSegment[US_STOP]) ? false : true;
		}
	}

	if (alignToSegment[MID]){		
		good_alignment = QuickAlignCheck(MID);//v1
		//good_alignment = false; //v2
		if (good_alignment){
			problemChild = false;
		}
		else{
			foundHit = FindDiagInSpecificRegion(MID);
			problemChild = (!foundHit) && (repeatEnd || repeatFwd) ? true : false;
		}
	}

	if (problemChild){
		minDiag = algnParams.fftParams.maxGap;
		maxDiag = algnParams.fftParams.maxGap;
	}
	else{
		if (repeatFwd){
			good_alignment = QuickAlignCheck(BEG_DS);//v1
			//good_alignment = false;//v2
			if (!good_alignment){
				foundHit = FindDiagInSpecificRegion(BEG_DS);
				if (!foundHit){
					minDiag = algnParams.fftParams.maxGap;
					maxDiag = algnParams.fftParams.maxGap;
				}
			}
		}

		if (repeatEnd){
			good_alignment = QuickAlignCheck(US_STOP);//v1
			//good_alignment = false;//v2
			if (!good_alignment){
				foundHit = FindDiagInSpecificRegion(US_STOP);
				if (!foundHit){
					minDiag = algnParams.fftParams.maxGap;
					maxDiag = algnParams.fftParams.maxGap;
				}
			}
		}
	}
	if (minDiag + maxDiag == 0){
		containsGaps = false;
		GaplessAlignment(); //no gaps in this cluster sequence, so align it along the diagonal
	}
	else{
		//SWAlignAndTrace();
		OverlapAlignAndTrace();
		containsGaps = true;
	}
}


void GermlineCluster::SearchGapsInFrontandEnd(){
	numSeqIns[0] = 0;
	numSeqIns[1] = 0;

	numSeqDel[0] = 0;
	numSeqDel[1] = 0;

	bool foundHit=false, foundHitEnd=false,  problemChild = false;
	querySeqToBeModified = (*queryInfoPointer).seq;

	DivideSeqIntoPeptides();
	
	foundHit = QuickAlignCheck(BEG);	//v1
	//foundHit = false;//v2
	if (!foundHit)
		foundHit = FindDiagInSpecificRegion(BEG);
		
	if (alignToSegment[LAST]){		
		foundHitEnd = QuickAlignCheck(LAST);//v1
		//foundHitEnd = false;//v2
		if (!foundHitEnd)
			foundHitEnd = FindDiagInSpecificRegion(LAST);			
	}
	
	problemChild = (!foundHit || !foundHitEnd);

	if (problemChild){
		minDiag = algnParams.fftParams.maxGap;
		maxDiag = algnParams.fftParams.maxGap;
	}
	
	if (minDiag + maxDiag == 0){
		containsGaps = false;
		GaplessAlignment(); //no gaps in this cluster sequence, so align it along the diagonal
	}
	else{
		//SWAlignAndTrace();
		OverlapAlignAndTrace();
		containsGaps = true;
	}
}



void GermlineCluster::GaplessAlignment(){
	/*int alignCoordinates[3];
	alignCoordinates[0] = alignmentToSequence[1];
	alignCoordinates[1] = alignmentToSequence[3];
	alignCoordinates[2] = alignmentToSequence[5];
	int SWScores[2];
	SWScores[0] = algnParams.swParams.matchScore;
	SWScores[1] = algnParams.swParams.mismatchScore;
	SWScore = dnafunctions::AlignDiagonal(alignCoordinates, querySeqToBeModified, paddedNTSequence, SWScores);
	finalScore = SWScore;
	*/
	int alignCoordinates[3];
	alignCoordinates[0] = alignmentToSequence[1];
	alignCoordinates[1] = alignmentToSequence[3];
	alignCoordinates[2] = alignmentToSequence[5];
	double gaplessScore;
	(*sw_align_cluster).GaplessAlignDiag(querySeqToBeModified, paddedNTSequence, alignCoordinates);
	gaplessScore = (*sw_align_cluster).maxScore;
	alignmentToSequence[1] = (*sw_align_cluster).seqStart;
	alignmentToSequence[2] = (*sw_align_cluster).seqEnd;
	alignmentToSequence[3] = (*sw_align_cluster).germStart;
	alignmentToSequence[4] = (*sw_align_cluster).germEnd;
	alignmentToSequence[5] = alignmentToSequence[2] - alignmentToSequence[1] + 1;
	finalScore = gaplessScore;
	SWgapPenalty = 0;
}

void GermlineCluster::GaplessAlignmentofGappedSeq(){
	int alignCoordinates[3];
	containsGaps = true;
	alignCoordinates[0] = alignmentToSequence[1];
	alignCoordinates[1] = alignmentToSequence[3];
	alignCoordinates[2] = alignmentToSequence[5];
	double gaplessScore;
	/*int SWScores[2];	
	SWScores[0] = algnParams.swParams.matchScore;
	SWScores[1] = algnParams.swParams.mismatchScore;
	gaplessScore = dnafunctions::AlignDiagonal(alignCoordinates, querySeqToBeModified, paddedNTSequence, SWScores);*/
	(*sw_align_cluster).GaplessAlignDiag(querySeqToBeModified, paddedNTSequence, alignCoordinates);
	gaplessScore = (*sw_align_cluster).maxScore;
	alignmentToSequence[1] = (*sw_align_cluster).seqStart;
	alignmentToSequence[2] = (*sw_align_cluster).seqEnd;
	alignmentToSequence[3] = (*sw_align_cluster).germStart;
	alignmentToSequence[4] = (*sw_align_cluster).germEnd;
	alignmentToSequence[5] = alignmentToSequence[2] - alignmentToSequence[1] + 1;
	SWScore = gaplessScore - SWgapPenalty;
	finalScore = SWScore;
}

void GermlineCluster::SWAlignAndTrace(){
	int SWcoordinates[7];
	string query = (*queryInfoPointer).seq;
	SWcoordinates[0] = alignmentToSequence[3];
	SWcoordinates[1] = alignmentToSequence[4];
	SWcoordinates[2] = alignmentToSequence[1] - alignmentToSequence[3];
	SWcoordinates[3] = minDiag;
	SWcoordinates[4] = maxDiag;
	SWcoordinates[5] = querySeqStart;
	SWcoordinates[6] = clusterSeqLength;
	(*sw_align_cluster).SWAlignDiag(query, paddedNTSequence, SWcoordinates);
	(*sw_align_cluster).Traceback();
	memcpy(numSeqIns, (*sw_align_cluster).numSeqIns, 2);
	memcpy(numSeqDel, (*sw_align_cluster).numSeqDel, 2);
	alignmentToSequence[1] = (*sw_align_cluster).seqStart;
	alignmentToSequence[2] = (*sw_align_cluster).seqEnd + numSeqIns[0] - numSeqDel[0]; //modify the end sequence to match the removed gaps and inserted insertion
	alignmentToSequence[3] = (*sw_align_cluster).germStart;
	alignmentToSequence[4] = (*sw_align_cluster).germEnd;
	alignmentToSequence[5] = alignmentToSequence[2] - alignmentToSequence[1]+1;
	querySeqToBeModified = (*sw_align_cluster).gaplessQuery;
	SWgapPenalty = (*sw_align_cluster).gapPenalties;
	SWScore = (*sw_align_cluster).maxScore;
	finalScore = SWScore;
}

void GermlineCluster::OverlapAlignAndTrace(){
	int SWcoordinates[7];
	string query = (*queryInfoPointer).seq;
	SWcoordinates[0] = alignmentToSequence[3];
	SWcoordinates[1] = alignmentToSequence[4];
	SWcoordinates[2] = alignmentToSequence[1] - alignmentToSequence[3];
	SWcoordinates[3] = minDiag;
	SWcoordinates[4] = maxDiag;
	SWcoordinates[5] = querySeqStart;
	SWcoordinates[6] = clusterSeqLength;
	(*sw_align_cluster).OverlapAlignDiag(query, paddedNTSequence, SWcoordinates);
	(*sw_align_cluster).Traceback();
	memcpy(numSeqIns, (*sw_align_cluster).numSeqIns, 2);
	memcpy(numSeqDel, (*sw_align_cluster).numSeqDel, 2);
	alignmentToSequence[1] = (*sw_align_cluster).seqStart;
	alignmentToSequence[2] = (*sw_align_cluster).seqEnd + numSeqIns[0] - numSeqDel[0]; //modify the end sequence to match the removed gaps and inserted insertion
	alignmentToSequence[3] = (*sw_align_cluster).germStart;
	alignmentToSequence[4] = (*sw_align_cluster).germEnd;
	alignmentToSequence[5] = alignmentToSequence[2] - alignmentToSequence[1] + 1;
	querySeqToBeModified = (*sw_align_cluster).gaplessQuery;
	SWgapPenalty = (*sw_align_cluster).gapPenalties;
	SWScore = (*sw_align_cluster).maxScore;
	finalScore = SWScore;
}

void GermlineCluster::SWAlignAndTrace(const string & query){
	int SWcoordinates[7];
	int qs = 0;
	while (query[qs] == NOT_NT && qs < query.length())
		qs++;
	SWcoordinates[0] = alignmentToSequence[3];
	SWcoordinates[1] = alignmentToSequence[4];
	SWcoordinates[2] = alignmentToSequence[1];
	SWcoordinates[3] = minDiag;
	SWcoordinates[4] = maxDiag;
	SWcoordinates[5] = qs;
	SWcoordinates[6] = clusterSeqLength;
	(*sw_align_cluster).SWAlignDiag(query, paddedNTSequence, SWcoordinates); //run local alignment using suggested diagonals
	(*sw_align_cluster).Traceback(); //traceback the local alignment result
	memcpy(numSeqIns, (*sw_align_cluster).numSeqIns, 2);
	memcpy(numSeqDel, (*sw_align_cluster).numSeqDel, 2);
	alignmentToSequence[1] = (*sw_align_cluster).seqStart;
	alignmentToSequence[2] = (*sw_align_cluster).seqEnd + numSeqIns[0] - numSeqDel[0]; //modify the end sequence to match the removed gaps and inserted insertion
	alignmentToSequence[3] = (*sw_align_cluster).germStart;
	alignmentToSequence[4] = (*sw_align_cluster).germEnd;
	querySeqToBeModified = (*sw_align_cluster).gaplessQuery;
	SWgapPenalty = (*sw_align_cluster).gapPenalties;
	SWScore = (*sw_align_cluster).maxScore;
	finalScore = SWScore;
}

void GermlineCluster::OverlapAlignAndTrace(const string & query){
	int SWcoordinates[7],qs =0;
	while (query[qs] == NOT_NT && qs < query.length())
		qs++;
	SWcoordinates[0] = alignmentToSequence[3];
	SWcoordinates[1] = alignmentToSequence[4];
	SWcoordinates[2] = alignmentToSequence[1];
	SWcoordinates[3] = minDiag;
	SWcoordinates[4] = maxDiag;
	SWcoordinates[5] = qs;
	SWcoordinates[6] = clusterSeqLength;
	(*sw_align_cluster).OverlapAlignDiag(query, paddedNTSequence, SWcoordinates); //run local alignment using suggested diagonals
	(*sw_align_cluster).Traceback(); //traceback the local alignment result
	memcpy(numSeqIns, (*sw_align_cluster).numSeqIns, 2);
	memcpy(numSeqDel, (*sw_align_cluster).numSeqDel, 2);
	alignmentToSequence[1] = (*sw_align_cluster).seqStart;
	alignmentToSequence[2] = (*sw_align_cluster).seqEnd + numSeqIns[0] - numSeqDel[0]; //modify the end sequence to match the removed gaps and inserted insertion
	alignmentToSequence[3] = (*sw_align_cluster).germStart;
	alignmentToSequence[4] = (*sw_align_cluster).germEnd;
	querySeqToBeModified = (*sw_align_cluster).gaplessQuery;
	SWgapPenalty = (*sw_align_cluster).gapPenalties;
	SWScore = (*sw_align_cluster).maxScore;
	finalScore = SWScore;
}

void GermlineCluster::FindAllGapsUsingAnotherRefCluster(const GermlineCluster & refCluster){
	bool problemChild = false, good_alignment, good_alignment_end;
	int ref_cluster_id = refCluster.clusterIndex;
	numSeqIns[0] = 0;
	numSeqDel[0] = 0;
	numSeqIns[1] = 0;
	numSeqDel[1] = 0;
	if (!refCluster.containsGaps){ //this cluster has no predicted gaps with query as determined by FFT
		if (alignmentToOtherClusters[ref_cluster_id].totalMutations != 0)
			LocateInsDelFromRef(refCluster);
		if (!containsGaps && alignmentToOtherClusters[ref_cluster_id].fftMutations == 0){ //for clustesr that have no intitali gaps in the sequence, and that the alignment to the reference cluster also does not have gaps that would be found by fft search (basically there should be no gaps anywhere)
			minDiag = alignmentToOtherClusters[ref_cluster_id].minDiag;
			maxDiag = alignmentToOtherClusters[ref_cluster_id].maxDiag;
			if (numSeqIns[0] + numSeqDel[0] > 0){ //count number of gaps, make sure they are more than 0
				FixGaps();
				good_alignment = QuickAlignCheck(BEG, refCluster); //check alignment using diagonal of best alignment
				good_alignment_end = QuickAlignCheck(LAST, refCluster);
				if (good_alignment && good_alignment_end)
					memcpy(alignmentToSequence, refCluster.alignmentToSequence, 6 * sizeof(int));
				else if (!good_alignment && !good_alignment_end){
					good_alignment = QuickAlignCheck(BEG); //check alignment using diagonal of best alignment
					good_alignment_end = QuickAlignCheck(LAST);
					if (!good_alignment || !good_alignment_end)
						problemChild = true;
				}
				else{
					problemChild = true;
				}

				//if (good_alignment)
				//	memcpy(alignmentToSequence, refCluster.alignmentToSequence, 6 * sizeof(int));
				//else{
				//	good_alignment = QuickAlignCheck(BEG); //if it ddint work, then try using diagonal of original alignment/this cluster
				//	problemChild = good_alignment ? problemChild : true;
				//}

				//good_alignment = QuickAlignCheck(LAST);

				//problemChild = good_alignment ? problemChild : true;

				if (!problemChild){
					GaplessAlignmentofGappedSeq();
				}
				else{
					minDiag = algnParams.swParams.maxGap;
					maxDiag = algnParams.swParams.maxGap;
					//SWAlignAndTrace();
					OverlapAlignAndTrace();
				}
			}
			else{
				querySeqToBeModified = queryInfoPointer->seq;
				GaplessAlignment();
			}
		}
		else{ //we are alignign the query to a sequence where FFT predicted has gaps
			minDiag = max(minDiag, alignmentToOtherClusters[ref_cluster_id].minDiag);
			maxDiag = max(maxDiag, alignmentToOtherClusters[ref_cluster_id].maxDiag);

			if (alignmentToOtherClusters[ref_cluster_id].fftMutations == numFFTGaps) //this means the gaps we observed were expected with repsect to the best cluster
			{
				FixGaps();
				good_alignment = QuickAlignCheck(BEG, refCluster);
				if (!good_alignment){
					good_alignment = QuickAlignCheck(BEG_DS, refCluster);
				}
				good_alignment_end = QuickAlignCheck(LAST, refCluster);
				if (!good_alignment_end)
					good_alignment_end = QuickAlignCheck(US_STOP, refCluster);
				if (!good_alignment_end && !good_alignment){
					good_alignment = QuickAlignCheck(BEG);
					if (!good_alignment)
						good_alignment = QuickAlignCheck(BEG_DS);
					good_alignment_end = QuickAlignCheck(LAST);
					if (!good_alignment_end)
						good_alignment_end = QuickAlignCheck(US_STOP);
					if (!good_alignment || !good_alignment_end){
						minDiag += 2;
						maxDiag += 2;
						//SWAlignAndTrace();						
						OverlapAlignAndTrace();
					}
					else{
						GaplessAlignmentofGappedSeq();
					}
				}
				else{
					memcpy(alignmentToSequence, refCluster.alignmentToSequence, 6 * sizeof(int));
					GaplessAlignmentofGappedSeq();
				}
			}
			else{
				bool repeatFwd, foundHit, repeatEnd;
				DivideSeqIntoPeptides();
				foundHit = FindDiagInSpecificRegion(BEG);
				repeatFwd = foundHit ? false : true;
				foundHit = FindDiagInSpecificRegion(LAST);
				repeatEnd = foundHit ? false : true;
				if (repeatFwd){
					foundHit = FindDiagInSpecificRegion(BEG_DS);
					if (!foundHit){
						minDiag = max(minDiag, algnParams.swParams.maxGap);
						maxDiag = max(maxDiag, algnParams.swParams.maxGap);
					}
				}
				if (repeatEnd){
					foundHit = FindDiagInSpecificRegion(US_STOP);
					if (!foundHit){
						minDiag = max(minDiag, algnParams.swParams.maxGap);
						maxDiag = max(maxDiag, algnParams.swParams.maxGap);
					}
				}

				//SWAlignAndTrace();
				OverlapAlignAndTrace();
			}
		}
	}
	else{
		if (minDiag + maxDiag == alignmentToOtherClusters[ref_cluster_id].minDiag + alignmentToOtherClusters[ref_cluster_id].maxDiag){ //this means the alignment is as expected. that is the predited FFT diags (mindiag and max diag) match the ins/deletion events predicted by aligning this cluster to the reference cluster
			minDiag = minDiag + refCluster.minDiag;//refCluster.minDiag -> these are ins/del events that we found between the equence and the BEST cluster. so we need to add in these posssible
			maxDiag = maxDiag + refCluster.maxDiag;
		}
		else{ //this means the min or max diags predicted from pairwise alignments do not match precisely. something else may be going on, so try to choose the maximum diagonal distance to align
			minDiag = refCluster.minDiag + max(minDiag,max(alignmentToOtherClusters[ref_cluster_id].minDiag, alignmentToOtherClusters[ref_cluster_id].maxDiag));
			maxDiag = refCluster.maxDiag + max(maxDiag,max(alignmentToOtherClusters[ref_cluster_id].minDiag, alignmentToOtherClusters[ref_cluster_id].maxDiag));				
		}
		//memcpy(alignmentToSequence, refCluster.alignmentToSequence, 6 * sizeof(int));
		//SWAlignAndTrace();
		OverlapAlignAndTrace();
	}
}

void GermlineCluster::FixGaps(){ //remove insertions and deletions in selected region of sequence
	int insPoints = 0, delPoints = 0;
	querySeqToBeModified = (*queryInfoPointer).seq;

	for (int j = 0; j < numSeqIns[1]; j++)
	{
		querySeqToBeModified.insert(seqIns[j][0] + insPoints + 1, seqIns[j][1], NOT_NT); //insert command works by inserting to the LEFT of the  position provided. our insertions are defeined as "position right before gap"
		insPoints += seqIns[j][1];
		for (int z = 0; z<numSeqDel[1]; z++){
			if (seqDel[z][0]>seqIns[j][0])
				seqDel[z][0] += insPoints;
		}
	}

	for (int i = 0; i < numSeqDel[1]; i++){
		querySeqToBeModified.erase(seqDel[i][0] - delPoints, seqDel[i][1]);
		delPoints += seqDel[i][1];
	}
}

void GermlineCluster::LocateInsDelFromRef(const GermlineCluster & refCluster){
	int temp, ref_cluster_id = refCluster.clusterIndex;
	SWgapPenalty = 0;
	numSeqIns[0] = 0;
	numSeqDel[0] = 0;
	numSeqIns[1] = 0;
	numSeqDel[1] = 0;
	for (int k = 0; k < alignmentToOtherClusters[ref_cluster_id].insertionEvents.size(); k++){
		temp = alignmentToOtherClusters[ref_cluster_id].insertionEvents[k][0] + (refCluster.alignmentToSequence[1] - refCluster.alignmentToSequence[3]);
		if (temp >  refCluster.alignmentToSequence[1] && temp < refCluster.alignmentToSequence[2]){ //important: make sure that the insertions for the pairwise clusters are actually present in the sequence /query of interest(i.e. if query is truncated then no need to correct mutations)
			numSeqIns[1]++;
			seqIns[numSeqIns[1] - 1][0] = temp;
			seqIns[numSeqIns[1] - 1][1] = alignmentToOtherClusters[ref_cluster_id].insertionEvents[k][1];
			numSeqIns[0] += seqIns[numSeqIns[1] - 1][1];
			SWgapPenalty += algnParams.swParams.swGapOpen + algnParams.swParams.swExtendGap*(seqIns[numSeqIns[1] - 1][1] - 1);
		}
	}

	int prevDelPos = -2;
	for (int k = 0; k < alignmentToOtherClusters[ref_cluster_id].deletionEvents.size(); k++){
		temp = alignmentToOtherClusters[ref_cluster_id].deletionEvents[k][0] + (refCluster.alignmentToSequence[1] - refCluster.alignmentToSequence[3]);
		//for (int l = 0; l < other_aligned_cluster_info->deletionEvents[k]; l++){
		if (temp>refCluster.alignmentToSequence[1] && temp < refCluster.alignmentToSequence[2]){ //important: make sure that the insertions for the pairwise clusters are actually present in the sequence /query of interest(i.e. if query is truncated then no need to correct mutations)
			if (temp - prevDelPos == 1){
				seqDel[numSeqDel[1] - 1][1] += 1;
				SWgapPenalty += algnParams.swParams.swExtendGap;
			}
			else{
				numSeqDel[1]++;
				seqDel[numSeqDel[1] - 1][0] = temp;// other_aligned_cluster_info->deletionEvents[k][0];
				seqDel[numSeqDel[1] - 1][1] = 1;// other_aligned_cluster_info->deletionEvents[k][1];
				SWgapPenalty += algnParams.swParams.swGapOpen;
			}
			prevDelPos = temp;// other_aligned_cluster_info->deletionEvents[k][0];
			numSeqDel[0] ++;
		}
	}
}