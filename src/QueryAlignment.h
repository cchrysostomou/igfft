#ifndef QUERYALIGN_H
#define QUERYALIGN_H

#include <algorithm>
#include <vector>
#include "StructVars.h"
#include "ReadFiles.h"
#include "DnaFunctions.h"
#include "GermlineCluster.h"
#include "FFTAlign.h"
#include "SWAlignment.h"
#include "string.h"
#include <fstream> //for ifstream

class QueryAlignment{
public:
	//QueryAlignment();
	QueryAlignment(const structvars::Fileinfo &, const structvars::AlignmentProgramSettings &);
	QueryAlignment(std::string &, const std::string &, const structvars::AlignmentProgramSettings &);
	QueryAlignment(const structvars::Fileinfo &, const std::string &, const std::string &, const structvars::AlignmentProgramSettings &);
	QueryAlignment(std::string &, const std::string &, const std::string &, const std::string &, const structvars::AlignmentProgramSettings &);
	~QueryAlignment();

	structvars::Fileinfo ReturnQueryFileInfo();
	structvars::Fileinfo GermlineQueryFileInfo();
	structvars::Fileinfo ClusterQueryFileInfo();

	void ResetQueryInfo(const structvars::FASTQinfo &, int, int);//for reading in new sequences/queries.  update the variable
	void ResetQueryInfo(structvars::ABRead *);//for reading in new sequences/queries.  update the variable

	//overloaded functions for different ways for reading in a germline database
	void ReadGermlineDatabase(const std::vector<std::string> &, int, bool allowSubstrings=true);		
	void ReadGermlineDatabase(const std::string &, const std::string &, bool allowSubstrings = true); //NO LONGER A FUNCTIONAL FUNCTION//
	int GetMinimumAlignmentLength();
	int FindBestAlignmentDiagonal(int guessDir = 0, const std::string & locus="");
	void OptimizeVGeneClusterAlignmentWithGaps(int); //this is a more exact alignmetn method to try to idnetify the actruall regions of alignment //to be called AFTER FindBestAlignmentDiagonal has been used. it wont work otherwise
	void OptimizeJGeneClusterAlignmentWithGaps(); 
	int AlignClusterResultsToGermline();
	int NumberUniqueLoci();
	int NumberUniqueChain();
	int UniqueLociListSize();//return maximum possilbe number of unique loci
	int UniqueChainListSize();//resutrn maximum possilbe number of unique chains
	void AddUniqueLociToVectorString(std::vector<std::string> &, int & currentLen);//if you pass in a vector of strings, then it updates this vector with unique Loci
	std::string UniqueLociString();
	std::string UniqueChainString();
	std::string FirstLoci();
	std::string FirstChain();
	
	void AlignToVGermlineMethod();
	void AlignToJGermlineMethod(); //this is different from Align to a v germline because we use the general assumptions that J germlines are usually shorter sequences and smaller, so the algorithm may be simpler

	void ModelOverlap();//this is the default model, nothing is sent in so we manually define the model between sequence and germline
	void ModelOverlap(const std::vector<int> &, const std::vector<int> &); //overloaded function if user wants to chnage the model
	//void SummarizeVGeneResults(std::map<std::string, std::string > &);
	//void SummarizeJGeneResults(std::map<std::string, std::string > &);
	void SummarizeVGeneResults(structvars::VGeneResults &);
	void SummarizeJGeneResults(structvars::JGeneResults &);
	//void SummarizeJGeneResults(std::map<std::string, std::string > &);

	void PrintClusterResults(const std::string &);

	clock_t ReturnTime();
	clock_t ReturnTimeOp();

private:
	//bool Compare(const std::vector<double>&, const std::vector<double> &);
	void MapCrossCorrelationToShift();
	void ReadGermlineFile(const std::string &, std::map<std::string, int> &, std::map<std::string, int> &);
	void ReadInClusterDatabase(const std::string &);//read a cluster file (tab delimited) into the proper variable
	void GroupGermlinesIntoClusters();
	void FindConsensusSeqInCluster(std::vector<int> &, const std::vector<std::vector<double>> &, std::vector<std::vector<int>>&, std::vector<std::string>&);
	void AlignPairwiseClusters(); //once the clusters have been craeted, we have to align each cluster to one another to learn how they align together and wehther difference clusters containg different gaps and where
	void FindGermlineStart(); //we also need to see if the germline that is within a cluster starts at the same position. if not, go ahead and pad the germline sequence with "N" 
	void AlignGermlineToQuery();
	void AnnotateQuery(int, const std::string &);

	void CopyGermlineDBToClusterDB();
	void InitializeTransformVariables(bool);
	void InitializeQueryVars();
	void InitializeSeqPeptideVar();
	void CopyQueryPrivateVarsToGermlineClusterVars();

	void CreateFFTSequenceDatabase();
	void DeleteGermlineVars();
	void DeleteClusterVars();
	void DeleteSeqPeptideVars();
	
	void DebugAlignments(int);//comment me out
	void UpdateUniqueLocus(const std::string & locus);//summarize the nubmer of unique loci found in results
	void UpdateUniqueChain(const std::string & locus);//summarize the number of unique chains found in results

	FFTAlign *full_seq_alignment, *full_seq_alignment_RC, *peptide_seq_alignment, ***list_of_fft_plans;
	std::vector<int> fftAlgnLens;
	int maxGermlineLength, maxQueryLength,maxClusterLength,numFFTPlans,currentPlan;
	
	structvars::ABRead querySeq[2];

	std::vector<std::string> annotationFields;

	int numUniqueChainHits, numUniqueLocusHits;
	std::vector<std::string> uniqueChains, uniqueLocus;
	int queryStart, queryEnd,germStart, germEnd, queryLen,seqShift;
	int numClusters, numGermlines, numClusterHits, numGermlineHits, num_with_gaps, num_without_gaps;
	structvars::Fileinfo germlineSeqInfo, clusterSeqInfo, queryInfo;
	std::vector<structvars::Germline> germlineDatabase;
	structvars::GermlineAlignmentResults *germlineHits;
	structvars::AlignmentProgramSettings algnSettings;
	complex_t **seqPeptideList; //will store all sequence peptides for FFT of small peptides
	int *seqPosStore; //will store indexes of what sequences have been currently transformed to FFT and what position need to be transformed
	double maxCorrectedScore,cluster_cutoff, minFFTAlgnScore;
	std::ofstream debugfile;
	GermlineCluster *aligned_to_clusters, **results_with_gaps, **results_no_gaps, *best_cluster_so_far;
	SWAlignment *finalized_sw_alignments;
	std::string *germlineSequences;
	std::vector<std::vector<double>> aligned_to_germline_scores;
	char bestSeqDir[2];
	clock_t s1, e1, t,tO;
};

inline void QueryAlignment::UpdateUniqueLocus(const std::string & locus ){
	bool newLocus = true;	
	for (int i = 0; i < numUniqueLocusHits;i++){
		if (locus == uniqueLocus[i]){
			newLocus = false;
			break;
		}
	}
	if (newLocus && locus!=""){
		uniqueLocus[numUniqueLocusHits] = locus;
		numUniqueLocusHits++;
	}	
	//if (numUniqueLocusHits>1)
		//printf("k");
}
inline void QueryAlignment::UpdateUniqueChain(const std::string & chain){
	bool newChain = true;	
	for (int i = 0; i < numUniqueChainHits;i++){
		if (chain == uniqueChains[i]){
			newChain = false;
			break;
		}
	}
	if (newChain && chain!=""){
		uniqueChains[numUniqueChainHits] = chain;
		numUniqueChainHits++;
	}
}

inline int QueryAlignment::NumberUniqueLoci(){
	return numUniqueLocusHits;
}
inline int QueryAlignment::NumberUniqueChain(){
	return numUniqueChainHits;
}
inline std::string QueryAlignment::UniqueChainString(){
	std::string chain = "";
	for (int i = 0; i < numUniqueChainHits - 1; i++)
		chain += uniqueChains[i] + ",";
	chain += uniqueChains[numUniqueChainHits - 1];
	return chain;
}

inline std::string QueryAlignment::UniqueLociString(){
	std::string loci = "";
	for (int i = 0; i < numUniqueLocusHits - 1; i++)
		loci += uniqueLocus[i] + ",";
	loci += uniqueLocus[numUniqueLocusHits - 1];
	return loci;
}

inline std::string QueryAlignment::FirstLoci(){
	return uniqueLocus[0];
}
inline std::string QueryAlignment::FirstChain(){
	return uniqueChains[0];	
}

inline int QueryAlignment::UniqueLociListSize(){//return maximum possilbe number of unique loci	
	return uniqueLocus.size();
}
inline int QueryAlignment::UniqueChainListSize(){ //resutrn maximum possilbe number of unique chains
	return uniqueChains.size();
}
inline int QueryAlignment::GetMinimumAlignmentLength(){ //resutrn maximum possilbe number of unique chains
	return minFFTAlgnScore;
}

#endif