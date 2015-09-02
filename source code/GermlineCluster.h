#ifndef GERMLINE_CLUSTER_H
#define GERMLINE_CLUSTER_H

#include <iomanip>
#include <iostream>
#include "string.h"

#include "FFTAlign.h"
#include "DnaFunctions.h"
#include "StructVars.h"
#include "SWAlignment.h"

class GermlineCluster{
	friend class QueryAlignment;
public:
	GermlineCluster();
	GermlineCluster(int, const std::string &, const structvars::AlignmentProgramSettings &);
	GermlineCluster(int, const std::string &, const structvars::AlignmentProgramSettings &, const std::vector<int> &, const std::vector<structvars::Germline> &, const std::vector<structvars::PairwiseClusterAlignment> &);
	~GermlineCluster();
	void InitializeComplex(const std::vector<int> &); //initalize the complex_t variables (allocate memory)
	void UpdateAlignmentToSequence(int *); //(inline function) modify the nucleotide alignments between the sequence and cluster by accounting for the true area where the sequence starts and cluster ends
	void FindAllGapsInQueryToClusterAlgn();
	void SearchGapsInFrontandEnd();
	void FindAllGapsUsingAnotherRefCluster(const GermlineCluster &);
	void GaplessAlignment();
	void GaplessAlignmentofGappedSeq();
	void SWAlignAndTrace();
	void SWAlignAndTrace(const std::string &); //overloaded fucntion, manually define the sequence
	void OverlapAlignAndTrace(); //perform overlap alignment instead of smithwaterman alignment
	void OverlapAlignAndTrace(const std::string &); //overloaded fucntion, manually define the sequence
	void LocateInsDelFromRef(const GermlineCluster &);
	void FixGaps();

	void DivideSeqIntoPeptides(); //just identifies coordinates to use for differenrt region within a peptide (top, bottom, middle, etc)
	bool FindDiagInSpecificRegion(region); //Will perform a FFT serach on a small peptide in the sequence designed by the region defined/given
	bool QuickAlignCheck(region); //will just count the number of matches between bases iwthin aspecific region of query
	bool QuickAlignCheck(region, const GermlineCluster &);

private:
	void AddGermlineMemberInfo(const std::vector<int> &, const std::vector<structvars::Germline> &);
	void AddPairwiseClusterInfo(const std::vector<structvars::PairwiseClusterAlignment> &);
	void AddSWAlignmentInfo(int, int);
	void MakePeptideFFT(int);

	//private variables
	FFTAlign *pointer_to_full_seq_alignment, *pointer_to_peptide_alignment;
	int **pointer_to_seq_peptide_list;
	complex_t ***pointer_to_fft_seq_peptides;
	std::string locusName;
	int clusterIndex, clusterSeqLength, numPeptides, numFFTPlans;
	structvars::AlignmentProgramSettings algnParams;
	structvars::ABRead *queryInfoPointer; //a pointer to the correct sequence /strand . contains nt sequence and complex sequence
	std::vector<structvars::PairwiseClusterAlignment> alignmentToOtherClusters;//this is an important variable, it stores how this cluster aligns to another cluster.  It stores the number of mutations and positions of insertion and deletions
	std::string consensus, paddedNTSequence, querySeqToBeModified;
	complex_t *complexSeqPad, **fourierSeq;
	complex_t **peptideFFT;

	int *germlineIndices, numGermlineMembers, *germlineStart;
	int algnPepResults[5][2][2]; //first dimension corresponds to each region (see enum region). the next dimension corresponds to either query or germilne (see enum seqtype). the last dimension corresponds to the starting nucleotide and to the ending nucleotide
	std::string *germlineSequences;
	SWAlignment *sw_align_cluster;
	bool alignToSegment[5];

	int **seqIns, **seqDel, numSeqIns[2], numSeqDel[2];

	int minDiag, maxDiag, alignmentToSequence[6], ds_us_distance,initialSWDiag;
	bool containsGaps;
	seqdir dir;
	int querySeqStart, querySeqEnd;
	int numFFTGaps, bestDiag;
	double FFTgappedScore, SWScore, finalScore, FFTgaplessScore, SWgapPenalty, minPerIdenPetideAlign;
	int numSteps, lineExit;
};

inline void GermlineCluster::UpdateAlignmentToSequence(int *seqAlignmentInfo){
	//structure of alignmentToSequence:position 0 => shift, 1=> nucleotide start position of sequence, 2=>nucleotide end position of sequence, 3=> nucleotide start position of germline, 4=> nucleotide end position of germline, 5=>alignment length
	memcpy(alignmentToSequence, seqAlignmentInfo, 6 * sizeof(int));
	for (int i = 1; i <= 4; i++)
		alignmentToSequence[i]--; //move the index of sequence/character position in string down by one to match index[0] = position one
	int shift = querySeqStart > alignmentToSequence[1] ? querySeqStart - alignmentToSequence[1] : 0; //if the starting position of the sequence is at a position higher than where the alignment to the querys equence begins, then shift the alignment to start where they strting squence is
	alignmentToSequence[1] += shift;
	alignmentToSequence[3] += shift;
	shift = clusterSeqLength-1 < alignmentToSequence[4] ? clusterSeqLength-1 - alignmentToSequence[4] : 0; //if the alignment position is a position greater than the end of the germline, then readjust position
	alignmentToSequence[2] += shift;
	alignmentToSequence[4] += shift;
	alignmentToSequence[5] = alignmentToSequence[2] - alignmentToSequence[1] + 1;
	initialSWDiag = alignmentToSequence[1] - alignmentToSequence[3];//best alignment diagonal predicted by fft
}

#endif
