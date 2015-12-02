#define _CRT_SECURE_NO_DEPRECATE

#ifndef _struct_vars_
#define _struct_vars_

#include <iostream>
#include <algorithm>
#include <fstream> //for ifstream
#include <cstring> //for strlen
#include <cstdio>
#include <vector>
#include <typeinfo> //for typeid
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string>
#include <map>
#include <sstream> //to use istringstream and getline
#include "fft.h"
//THIS FILE DEFINES ALL OF THE STRUCTURE VARIABLES WE USE IN THE PROGRAM

enum alignment_movement { DIAG, RIGHT, DOWN, START }; //DEFINES THE ORDER IN WHICH WE PREFER THE SW ALIGNMENT PATH (I.E. PREFER MOVING DIAGONALLY BEFORE INTRODUCING INSERTIONS/DELETIONS)
enum direction{ FRONT, END };
enum tracebackInfo{MATCHES,MISMATCHES,INDEL,SCORE,GAPPENALTY,SEQPOS,ALGNPOS};
enum method{ NONE, COPY, MAKE };

enum region{ BEG, BEG_DS, MID, US_STOP, LAST }; //different regions
enum seqtype{ QUERY, GERMLINE }; //corresponds to query
enum position{ FIRST, FINAL }; //start and stop position of an alignment
enum seqdir{ FS,RC }; //start and stop position of an alignment

const char NOT_NT = 'N';

namespace structvars{
	//STRUCTURE FOR FASTA DATA
	struct FASTAinfo{
		std::string seqHeader, seq;
	};

	//STRUCTOR FOR FASTQ DATA
	struct FASTQinfo{
		std::string seqHeader, seq, quality;
		FASTQinfo(){
			quality = "";
		}
	};

	//STRUCTURE FOR FILE INFO
	struct Fileinfo{
		int numSeqs, maxSeqLen, minSeqLen,avgSeqLen;
		std::vector<std::vector<int>> lengthDistribution;
	};

	//just to store all of the possible sub fields of an antiobidy (i.e. cdr1 and 2)
	struct Abregion{
		std::string sequence;
		int startpos, endpos,startpos_algnseq,endpos_algnseq;
	};

	//antibody read
	struct ABRead{
		std::string seq, header, seq_rev;// , quality;
		int seqLength, seqStart;
		complex_t* complSeq;
		complex_t* complSeqPad;
		complex_t* fourSeq;
		ABRead(){
			complSeq = NULL;
			complSeqPad = NULL;
			fourSeq = NULL;
			seqStart = 0;
			//quality = "";
		}
		~ABRead(){
			if (complSeq != NULL)
				free(complSeq);			
			if (complSeqPad != NULL)
				free(complSeqPad);
			if (fourSeq != NULL)
				free(fourSeq);
			complSeq = NULL;
			complSeqPad = NULL;
			fourSeq = NULL;
		}
	};

	//defines teh structure of each type of germline sequence
	struct Germline{
		int id, seqLength, seqStart;
		ABRead germlineSeqData; //all the information we will need for aligning the sequences
		std::string genename, locus,chain;
		std::map<std::string, std::string> additional_info;
		std::map<std::string, Abregion> annotation; //dictionary for labeling a germline by 'FR1','FR2','CDR1','CDR2','CDR3', etc
		std::vector<int> annotationIndex; //each value in this array corresponds to the annotated region at each position. 

		Germline(){
			id = -1;
			seqLength = 0;
			seqStart = 0;
		}
	};

	//FFT settings
	struct FFTSettings{
		double fftGapOpen, fftExtendGap;
		int maxGap, extraGap;
		double scoreCutoff, sensitivity;
		double cluster_threshold;
		int guessDir;

		FFTSettings(){
			guessDir = 0;
			fftGapOpen = 10;
			fftExtendGap = 2;
			maxGap = 10;
			extraGap = 5;
			sensitivity = 5.0;
			scoreCutoff = 0.3;// 0.5;
			cluster_threshold = 0.8;
		};
	};

	//smith waterman settings
	struct SWAlignSettings{
		int matchScore; //this is for local/SW alignments
		int mismatchScore; //this is for local/SW alignments
		double swGapOpen, swExtendGap;
		int maxGap;
		SWAlignSettings(){
			matchScore = 5;
			mismatchScore = 4;
			swGapOpen = 50;
			swExtendGap = 10;
			maxGap = 10;
		}
	};

	struct PeptideFFTSettings{
		int peptide_len, length_flank_gap;
		double peptide_gap_sensitivity;
		PeptideFFTSettings(){
			peptide_len = 20;
			length_flank_gap = 5;
			peptide_gap_sensitivity = 3;
		}
	};

	struct AlignmentProgramSettings{
		std::vector<std::string> fileGermline; //current format for each string from user "filetype(fasta/txt),filename referring to location of germline sequences, filename referring to location of clustered sequences"
		FFTSettings fftParams;
		SWAlignSettings swParams;
		PeptideFFTSettings peptideParams;
		bool group_into_clusters;
		
		int maxHits, numObsAboveFoldRatio;
		double scoreFoldRatio, clusterGermlineCutoff;	
		AlignmentProgramSettings(){			
			group_into_clusters = true;
			scoreFoldRatio = 2.0;
			clusterGermlineCutoff = 0.8;
			maxHits = 5;
			numObsAboveFoldRatio = 3;		
		}
	};

	//data for alignment to other clusters (gap information)
	struct PairwiseClusterAlignment{ //information about an alignment between two sequences
		int minDiag, maxDiag, bestDiag;//stores minimm diagonal and maximum diagonal betwene an alignment
		int fftMutations, totalMutations;//stores mutations predicted by FFT and total mutations ???
		//int numAlgn;
		//int top_start; //start position of the cluster aligned??
		//int bottom_start; //start nucleotide position of the cluster it is aligned to???
		std::vector<std::vector<int>> insertionEvents; //stores the postion of each insertion and the number of times a consecutive insertion occurs
		std::vector<std::vector<int>> deletionEvents; //stores teh position of each deletion
	};
	
	struct GermlineAlignmentResults{
		double score;
		char dir;
		int index, numMismatch, numMatch, numIns, numDel,minDiag,maxDiag, clusterIndex;
		std::string germlineAlSeq, queryAlSeq;
		int startQuery, startGerm, endQuery, endGerm, algnLen, swAlignDiag;
				
		//int match[6], mismatch[6],indel[6];
		std::map<std::string, int> match;
		std::map<std::string, int> mismatch;
		std::map<std::string, int> indel;
		std::map<std::string, structvars::Abregion> annotation;
		std::map<std::string, bool> included;
		//structvars::Abregion annotation[6];		
		int algnRegion[2];//structure will be FR1/FR3 or CDR1/FR3, or FR1/CDR3 etc,etc
	};

	struct VGeneResults{
		char direction;
		int QS, QE, GS, GE, totalMatch, totalMismatch,totalIndel,cdr3start,germlineCodingFrame,codonStart;
		std::string startingAnnotation,endingAnnotation;
		std::string algnRegion[2];//structure will be FR1...FR3 or CDR1...FR3, or FR1...CDR3 etc,etc
		std::map<std::string, std::string> numMatch;
		std::map<std::string, std::string> numMisMatch;
		std::map<std::string, std::string> numIndel;
		std::map<std::string, std::string> sequence;
		std::map<std::string, std::string> startPosText; 
		std::map<std::string, int> startPos;
		std::map<std::string, std::string> endPos;

		std::map<std::string, std::string> startPos_algnseq;
		std::map<std::string, std::string> endPos_algnseq;


		std::map<std::string, std::string> readingFrameText;
		std::map<std::string, int> readingFrame;
		std::string genes,scores;	
		std::string queryAlSeq, germlineAlSeq;
	};

	struct JGeneResults{
		char direction;
		int QS, QE, GS, GE, totalMatch, totalMismatch, totalIndel,germlineCodingFrame,codonStart;
		std::string genes, scores;
		std::string queryAlSeq, germlineAlSeq;

	};


	//defines the input parameters for the program
	struct VarParameters{
		//define variables to use for program
		std::string fileinputName;
		std::string fileoutputName;
		std::string inputFileFormat;
		bool containsVGermline, containsJGermline, containsDGermline;
		std::map<char, structvars::AlignmentProgramSettings> query_algn_settings;
		bool ignoreHeader;
		VarParameters()
		{			
			inputFileFormat = "FASTA";
			ignoreHeader = true;
			containsVGermline = false;
			containsJGermline = false;
			containsDGermline = false;
			query_algn_settings['V'] = AlignmentProgramSettings();
			query_algn_settings['D'] = AlignmentProgramSettings();
			query_algn_settings['J'] = AlignmentProgramSettings();
		};
	};

	

}
#endif