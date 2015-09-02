#ifndef SWALIGN_H
#define SWALIGN_H

#include "StructVars.h"
#include <iostream>
#include <iomanip>
#include <limits>

class SWAlignment{
	friend class QueryAlignment;
	friend class GermlineCluster;
public:
	SWAlignment(int, int);
	SWAlignment(int, int, structvars::SWAlignSettings);
	
	~SWAlignment();

	void printAlignment();
	void printTracebackAnnotation(int);
	void SWAlignDiag(const std::string &, const std::string &, int *);
	void OverlapAlignDiag(const std::string &, const std::string &, int *);
	void OverlapAlignComplete(const std::string &, const std::string &);
	void SWAlignComplete(const std::string &, const std::string &, int *);
	void GaplessAlignDiag(const std::string &, const std::string &, int *);
	void Traceback();
	void Traceback(int ***); //overloaded function: you can pass in another socring matrix and alignment path to perform traceback
	void TracebackWithAnnotation(); //this is a special traceback function used in the final alignment step. In addition to finding the traceback path, it will also store data for the number of matches, mismatches, indels, and query positions

private:
	structvars::SWAlignSettings swParams;
	int germLen, seqLen;
	double ***scoreMatrix, *gaplessScores,**tracebackData; //matrix of alignment scores, //for gapless alignments	
	int ***alignmentPath, final_posIJ[2];
	double maxScore, gapPenalties;
	int totalMatch, totalMismatch, totalIndel;
	int **seqDel, **seqIns, numSeqIns[2], numSeqDel[2]; //seqDel[row] = position of deletin, seqDel[column]=number of consecutive dleetions at that position //numSeqIns[0] => total insertions, numSeqDel[0] = > total deletions, numSeqIns[1] = > total insertion rows (number of times a new gap is opened), numSeqDel[1] => total deletion rows (number of times a new gap is opened)
	int seqStart, seqEnd, germStart, germEnd;
	int maxSeqInsSize, maxSeqDelSize;
	std::string query, germline, alignedGerm, alignedQuery, gaplessQuery, gaplessGermline;
};

#endif