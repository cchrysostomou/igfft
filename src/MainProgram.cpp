#define _CRT_SECURE_NO_DEPRECATE
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
#include <sstream> //to use istringstream and getline

#include "QueryAlignment.h"
#include "ReadFiles.h" //this file stores fucntions for reading files in an out
#include "StructVars.h" //this file stores all structure variables used in this program
#include "DnaFunctions.h" //for function invovling nucleotides
#include "FFTFunctions.h"
#include <ctime>
#include <limits>

using namespace std;

structvars::VarParameters variables_used; //stores useful variables for running program (all input variables by user get stored here)
enum UserSettings{
	allowed_gaps,
	match_sw,
	mismatch_sw,
	gap_open_sw,
	gap_extend_sw,
	gap_open_fft,
	gap_extend_fft,
	sensitivity_fft,
	sensitivity_peptide,
	similar_clusters,
	peptide_len,
	min_per_id_threshold,
	max_germline_hits,
	score_ratio,
	num_above_score_ratio,
	cluster_germlines_cutoff,
	group_clusters
};
struct settings_description{
	UserSettings parameter;
	string description;
	settings_description(){
		description = "";
	}
};
static map<string, settings_description> interpreter;

void InitializeUserSettings(){
	interpreter["gap"].parameter = allowed_gaps;
	interpreter["match_sw"].parameter = match_sw;
	interpreter["mismatch_sw"].parameter = mismatch_sw;
	interpreter["gap_open_sw"].parameter = gap_open_sw;
	interpreter["gap_extend_sw"].parameter = gap_extend_sw;
	interpreter["gap_open_fft"].parameter = gap_open_fft;
	interpreter["gap_extend_fft"].parameter = gap_extend_fft;
	interpreter["s_fft"].parameter = sensitivity_fft;
	interpreter["s_pep"].parameter = sensitivity_peptide;
	interpreter["similar_clusters"].parameter = similar_clusters;
	interpreter["pep_len"].parameter = peptide_len;
	interpreter["min_per_id"].parameter = min_per_id_threshold;
	interpreter["num_hits"].parameter = max_germline_hits;
	interpreter["ratio_cutoff"].parameter = score_ratio;
	interpreter["times_above_ratio"].parameter = num_above_score_ratio;
	interpreter["cluster_per_id"].parameter = cluster_germlines_cutoff;
	interpreter["group_clusters"].parameter = group_clusters;
}

map<string, int> output_file_fields;

void EvaluateParameters(int, char *[]); //This function will evaluate the parameters inputed from the user
void EvaluateAlgnParameters(char, string, string); //This function will evaluate the parameters inputed from the user
void RunAlignmentProgram(map<string, int>); //main function for aligning a set of sequences to the v germline
vector<string> WriteHeaderRow(ofstream &); //write the initial header row of the file
void WriteRow(ofstream & outfile, vector<string> fields, int); //write results to a row in the file
void SummarizeData(vector<string> &, const vector<structvars::Germline> & germline_set, const vector<vector<double>> & germline_scores, const int & numGenes); //organize results into a vector to be easily parsed by "writerow" function
double GetDecimalValue(string parameter, bool between_0_1);
int GetIntegerValue(string parameter, bool between_0_1);
void PrintDefaultParameters(bool to_file = true);

int main(int argc, char *argv[]) //read in variables and determined how we are aligning sequences to provided germlines
{
	InitializeUserSettings();
	map<string, int> alignmentMethod; //sturucture of alignmentMethodl -> {'vgene': '0, 1, or 2', 'jgene': '0, 1, or 2'} where 0 = > no alignment to germline, 1 => use clusters to align, 2 => use unclustered germlines to align
	alignmentMethod["vgene"] = NONE;
	alignmentMethod["jgene"] = NONE;
	alignmentMethod["dgene"] = NONE;

	//First, evaluate command line parameters inputed by the user	
	EvaluateParameters(argc, argv);

	//NEXT START IMPORTING ALL OF THE FILES AND READING IN THE GERMLINE DATABASE AND SEQUENCES //////
	//if a vgene database is passed in
	if (variables_used.containsVGermline){
		alignmentMethod["vgene"] = variables_used.query_algn_settings['V'].group_into_clusters ? MAKE : COPY;
	}

	//if a jgene database is passed in
	if (variables_used.containsJGermline){
		alignmentMethod["jgene"] = variables_used.query_algn_settings['J'].group_into_clusters ? MAKE : COPY;
	}

	//if a dgene database is passed in
	if (variables_used.containsDGermline){
		//REGARDLESS OF USER INPUT SETTING, ALWAYS COPY THE GERMLINE OF D GENES (DONT CLUSTER)
		alignmentMethod["dgene"] = COPY; //variables_used.query_algn_settings['D'].group_into_clusters ? MAKE : COPY;
	}

	RunAlignmentProgram(alignmentMethod);

	std::printf("\nAnalysis complete\n");
	return 0;
}

void RunAlignmentProgram(map<string, int> method)
{
	ofstream outfile(variables_used.fileoutputName); /*VARIABLES FOR WRITING THE RESULTS TO A FILE*/
	vector<string> headers = WriteHeaderRow(outfile);

	vector<string> results_to_text(headers.size());

	string notes;
	/*********************************************/

	structvars::FASTQinfo seqRead; //for reading in sequences (queries)
	structvars::Fileinfo inputseqinfo; //info on the sequence query file (i.e. num sequences, max seq length, etc)
	string seq_filename = variables_used.fileinputName; //name of the query sequence file

	//step one, read in the query file and determine: a) # of sequences and b) max sequence length
	if (variables_used.inputFileFormat == "FASTA"){
		std::printf("Calcuating number of sequences\n");
		inputseqinfo = readfiles::convertFASTAtoTAB(seq_filename);
	}
	else if (variables_used.inputFileFormat == "FASTQ"){
		std::printf("Calcuating number of sequences\n");
		inputseqinfo = readfiles::convertFASTQtoTAB(seq_filename);
	}
	else{
		std::printf("Calcuating number of sequences\n");
		inputseqinfo = readfiles::readSeqLengths(seq_filename.c_str(), "TAB", variables_used.ignoreHeader);
	}

	//Initialize a vgene alignment variable and a jgene alignment variable
	QueryAlignment vGermlineAlignment(inputseqinfo, variables_used.query_algn_settings['V']), jGermlineAlignment(inputseqinfo, variables_used.query_algn_settings['J']); //The two main variables used in the program //each of these variables will be used to align the queries to either the v gene database or j gene databse
	QueryAlignment dGermlineAlignment(inputseqinfo, variables_used.query_algn_settings['D']);

	if (method["vgene"] != NONE){ //First we need to read in the file containing v germline genes.  Reading in teh database willl subsequently merge germlinesinto clusters
		std::printf("Processing the provided V germline database:\n");
		if (variables_used.query_algn_settings['V'].fileGermline[1].length()>0){
			vGermlineAlignment.ReadGermlineDatabase(variables_used.query_algn_settings['V'].fileGermline[0], variables_used.query_algn_settings['V'].fileGermline[1], false); //read in the v germline sequences, and the second provided file defining the clustered sequences
		}
		else{
			vGermlineAlignment.ReadGermlineDatabase(variables_used.query_algn_settings['V'].fileGermline[0], method["vgene"], false); //no cluster was provided in a file, so we will cluster sequences ad-hoc
		}
		vGermlineAlignment.PrintClusterResults("vgermlineclusters.txt");
	}

	if (method["jgene"] != NONE){ //First we need to read in the file containing j germline genes.  Reading in teh database willl subsequently merge germlinesinto clusters
		std::printf("Processing the provided J germline database:\n");
		if (variables_used.query_algn_settings['J'].fileGermline[1].length() > 0){
			jGermlineAlignment.ReadGermlineDatabase(variables_used.query_algn_settings['J'].fileGermline[0], variables_used.query_algn_settings['J'].fileGermline[1]);
		}
		else{
			jGermlineAlignment.ReadGermlineDatabase(variables_used.query_algn_settings['J'].fileGermline[0], method["jgene"]);//no cluster was provided in a file, so we will cluster sequences ad-hoc
		}
		jGermlineAlignment.PrintClusterResults("jgermlineclusters.txt");
	}

	if (method["dgene"] != NONE){ //First we need to read in the file containing j germline genes.  Reading in teh database willl subsequently merge germlinesinto clusters
		std::printf("Processing the provided D germline database:\n");
		if (variables_used.query_algn_settings['D'].fileGermline[1].length() > 0){
			jGermlineAlignment.ReadGermlineDatabase(variables_used.query_algn_settings['D'].fileGermline[0], variables_used.query_algn_settings['D'].fileGermline[1]);
		}
		else{
			jGermlineAlignment.ReadGermlineDatabase(variables_used.query_algn_settings['D'].fileGermline[0], method["dgene"]);//no cluster was provided in a file, so we will cluster sequences ad-hoc
		}
	}

	/****************************************/
	int counts = 0, numVClusterHits = 0, numVGermlineHits = 0, numJClusterHits = 0, numJGermlineHits = 0;
	ifstream inFile(seq_filename.c_str()); //input filename containing teh query sequences
	int totalHits = 0;
	int startPer = 0, perIndiciate = 10, perDone = 0, seqAnalyzed = 0, total_seqs = inputseqinfo.numSeqs; //for reporting percent done/status of alignment
	double cluster_cutoff;
	int jStart, dStart, abstart, abend, cdr3start, v_start_temp, j_start_temp,codonStart,codonFrame;
	//std::map<std::string, std::string> VGeneResults;
	structvars::VGeneResults vGeneResults;
	structvars::JGeneResults jGeneResults;
	clock_t s1, s2, t = 0, tna = 0;
	bool vfound, jfound, dfound;
	char dir;
	string temp,locusSearch,uniqueLoci,uniqueChains;
	vector<string> uniqueLocusVector;
	uniqueLocusVector.resize(vGermlineAlignment.UniqueLociListSize() + jGermlineAlignment.UniqueLociListSize());
	int currentLocusNum; 
	if (variables_used.ignoreHeader && variables_used.inputFileFormat != "FASTA" && variables_used.inputFileFormat != "FASTQ")
		getline(inFile, temp); //skip the header line, if its not a fasta file

	clock_t total_time = clock(), current_clock;
	//clock_t ;
	double analysis_time;
	string startingAnnotation, endingAnnotation;
	printf("\nGermline Alignment Started\n");
	/*START OF THE ALIGNMENT PROGRAM*/
	int minVAlignment = vGermlineAlignment.GetMinimumAlignmentLength(), minJAlignment = jGermlineAlignment.GetMinimumAlignmentLength();
	//while (seqAnalyzed<inputseqinfo.numSeqs && !inFile.eof()){// !inFile.eof()){//seqAnalyzed < 10000){ //(seqAnalyzed<1000){// (!inFile.eof()){ //(seqAnalyzed < 5000){
	for(int seqAnalyzed=0;seqAnalyzed<inputseqinfo.numSeqs;seqAnalyzed++){
		if(inFile.eof())
			continue;		

		//cout<<seqAnalyzed<<endl;
		//cout<<inputseqinfo.numSeqs<<endl;		
		cluster_cutoff = 0;
		
		vfound = false;
		jfound = false;
		dfound = false;
		dir = '\0';

		/*print out the current status of the program/how many sequences have been analyzed */
		perDone = int((seqAnalyzed / float(total_seqs)) * 100);
		if (perDone%perIndiciate == 0 && perDone > startPer){
			current_clock = clock();
			analysis_time = ((current_clock - total_time)*1.0 / CLOCKS_PER_SEC) / 60;
			printf("%i Percent Complete: %g minutes \n", perDone, analysis_time);
			startPer = perDone;
		}
		/**/
		
		for (int k = 0; k < results_to_text.size(); k++)	
			results_to_text[k] = "";
		
		seqRead = readfiles::readIndividualTabSequence(inFile); //read in the query from the file
		
		if(seqRead.seqHeader=="Noneerror" && seqRead.seq=="")
			continue; //empty line provided
		else if (seqRead.seq == ""){
			results_to_text[output_file_fields["Header"]] = seqRead.seqHeader;
			results_to_text[output_file_fields["Sequence quality"]] = "";
			results_to_text[output_file_fields["Notes"]] = "Sequence not found";
			WriteRow(outfile, results_to_text, headers.size());			
			continue;
		}
		
		//seqAnalyzed++;
				
		notes = "";
		abstart = 0;
		abend = 0;
		vGeneResults.QS = 0;
		vGeneResults.QE = 0;
		codonStart = 0;
		
		vGeneResults.cdr3start = -1;
		startingAnnotation = "";
		endingAnnotation = "";
		currentLocusNum = 0;
		numVGermlineHits = 0;
		numJGermlineHits = 0;
		numVClusterHits = 0;
		numJClusterHits = 0;
		if (method["vgene"] != NONE && seqRead.seq.length()>minVAlignment){ //PERFORM ALIGNMENT TO VGENES
			vGermlineAlignment.ResetQueryInfo(seqRead, 0, seqRead.seq.length()); //now update relevant sequence variable we use for aligning with seqRead info
			s1 = clock();
			numVClusterHits = vGermlineAlignment.FindBestAlignmentDiagonal();
			s2 = clock();
			t += s2 - s1;
			if (numVClusterHits > 0){
				s1 = clock();
				vfound = true;
				vGermlineAlignment.OptimizeVGeneClusterAlignmentWithGaps(seqAnalyzed);
				s2 = clock();
				tna += s2 - s1;
				numVGermlineHits = vGermlineAlignment.AlignClusterResultsToGermline();

				if (numVGermlineHits > 0){
					//###summarize vgene resutls###/
					vGermlineAlignment.SummarizeVGeneResults(vGeneResults);
					//v_start_temp = vGeneResults.QS;

					//dnafunctions::AdjustToCorrectFrame(vGeneResults.QS, vGeneResults.GS, vGeneResults.germlineCodingFrame);
					abstart = vGeneResults.QS;
					abend = vGeneResults.QE;
					codonStart = vGeneResults.codonStart;
					dir = vGeneResults.direction;
					startingAnnotation = vGeneResults.startingAnnotation;
					endingAnnotation = vGeneResults.endingAnnotation;
					//vGeneResults.startPos[vGeneResults.startingAnnotation] = abstart;

					results_to_text[output_file_fields["Direction"]] = vGeneResults.direction;
					results_to_text[output_file_fields["Top_V-Gene_Hits"]] = vGeneResults.genes;
					results_to_text[output_file_fields["V-Gene_Alignment_Scores"]] = vGeneResults.scores;

					results_to_text[output_file_fields["VGENE: Query_Start"]] = std::to_string(abstart);
					results_to_text[output_file_fields["VGENE: Query_End"]] = std::to_string(vGeneResults.QE);
					results_to_text[output_file_fields["VGENE: Germline_Start"]] = std::to_string(vGeneResults.GS);
					results_to_text[output_file_fields["VGENE: Germline_End"]] = std::to_string(vGeneResults.GE);

					results_to_text[output_file_fields["FR1_Sequence.NT"]] = vGeneResults.sequence["FR1"];
					results_to_text[output_file_fields["FR2_Sequence.NT"]] = vGeneResults.sequence["FR2"];
					results_to_text[output_file_fields["FR3_Sequence.NT"]] = vGeneResults.sequence["FR3"];
					results_to_text[output_file_fields["CDR1_Sequence.NT"]] = vGeneResults.sequence["CDR1"];
					results_to_text[output_file_fields["CDR2_Sequence.NT"]] = vGeneResults.sequence["CDR2"];

					results_to_text[output_file_fields["FR1_Sequence.AA"]] = dnafunctions::TranslateSeq(vGeneResults.sequence["FR1"],vGeneResults.readingFrame["FR1"],vGeneResults.startPos["FR1"]);
					results_to_text[output_file_fields["FR2_Sequence.AA"]] = dnafunctions::TranslateSeq(vGeneResults.sequence["FR2"], vGeneResults.readingFrame["FR2"], vGeneResults.startPos["FR2"]);
					results_to_text[output_file_fields["FR3_Sequence.AA"]] = dnafunctions::TranslateSeq(vGeneResults.sequence["FR3"], vGeneResults.readingFrame["FR3"], vGeneResults.startPos["FR3"]);
					results_to_text[output_file_fields["CDR1_Sequence.AA"]] = dnafunctions::TranslateSeq(vGeneResults.sequence["CDR1"], vGeneResults.readingFrame["CDR1"], vGeneResults.startPos["CDR1"]);
					results_to_text[output_file_fields["CDR2_Sequence.AA"]] = dnafunctions::TranslateSeq(vGeneResults.sequence["CDR2"], vGeneResults.readingFrame["CDR2"], vGeneResults.startPos["CDR2"]);					

					results_to_text[output_file_fields["VGENE: Query_FR1_Start::End"]] = vGeneResults.startPosText["FR1"] + "::" + vGeneResults.endPos["FR1"];
					results_to_text[output_file_fields["VGENE: Query_FR2_Start::End"]] = vGeneResults.startPosText["FR2"] + "::" + vGeneResults.endPos["FR2"];
					results_to_text[output_file_fields["VGENE: Query_FR3_Start::End"]] = vGeneResults.startPosText["FR3"] + "::" + vGeneResults.endPos["FR3"];
					results_to_text[output_file_fields["VGENE: Query_CDR1_Start::End"]] = vGeneResults.startPosText["CDR1"] + "::" + vGeneResults.endPos["CDR1"];
					results_to_text[output_file_fields["VGENE: Query_CDR2_Start::End"]] = vGeneResults.startPosText["CDR2"] + "::" + vGeneResults.endPos["CDR2"];

					results_to_text[output_file_fields["VGENE: Alignment_FR1_Start::End"]] = vGeneResults.startPos_algnseq["FR1"] + "::" + vGeneResults.endPos_algnseq["FR1"];					
					results_to_text[output_file_fields["VGENE: Alignment_FR2_Start::End"]] = vGeneResults.startPos_algnseq["FR2"] + "::" + vGeneResults.endPos_algnseq["FR2"];
					results_to_text[output_file_fields["VGENE: Alignment_FR3_Start::End"]] = vGeneResults.startPos_algnseq["FR3"] + "::" + vGeneResults.endPos_algnseq["FR3"];
					results_to_text[output_file_fields["VGENE: Alignment_CDR1_Start::End"]] = vGeneResults.startPos_algnseq["CDR1"] + "::" + vGeneResults.endPos_algnseq["CDR1"];
					results_to_text[output_file_fields["VGENE: Alignment_CDR2_Start::End"]] = vGeneResults.startPos_algnseq["CDR2"] + "::" + vGeneResults.endPos_algnseq["CDR2"];

					results_to_text[output_file_fields["VGENE: Alignment_Sequence_Query"]] = vGeneResults.queryAlSeq;
					results_to_text[output_file_fields["VGENE: Alignment_Sequence_Germline"]] = vGeneResults.germlineAlSeq;

					results_to_text[output_file_fields["VGENE: Total_Matches"]] = std::to_string(vGeneResults.totalMatch);
					results_to_text[output_file_fields["VGENE: Total_Mismatches"]] = std::to_string(vGeneResults.totalMismatch);
					results_to_text[output_file_fields["VGENE: Total_Indel"]] = std::to_string(vGeneResults.totalIndel);

					results_to_text[output_file_fields["VGENE_Matches: FR1,CDR1,FR2,CDR2,FR3,CDR3"]] = vGeneResults.numMatch["FR1"] + "," + vGeneResults.numMatch["CDR1"] + "," + vGeneResults.numMatch["FR2"] + "," + vGeneResults.numMatch["CDR2"] + "," + vGeneResults.numMatch["FR3"] + "," + vGeneResults.numMatch["CDR3"];
					results_to_text[output_file_fields["VGENE_Mismatches: FR1,CDR1,FR2,CDR2,FR3,CDR3"]] = vGeneResults.numMisMatch["FR1"] + "," + vGeneResults.numMisMatch["CDR1"] + "," + vGeneResults.numMisMatch["FR2"] + "," + vGeneResults.numMisMatch["CDR2"] + "," + vGeneResults.numMisMatch["FR3"] + "," + vGeneResults.numMisMatch["CDR3"];
					results_to_text[output_file_fields["VGENE_Indels: FR1,CDR1,FR2,CDR2,FR3,CDR3"]] = vGeneResults.numIndel["FR1"] + "," + vGeneResults.numIndel["CDR1"] + "," + vGeneResults.numIndel["FR2"] + "," + vGeneResults.numIndel["CDR2"] + "," + vGeneResults.numIndel["FR3"] + "," + vGeneResults.numIndel["CDR3"];
					results_to_text[output_file_fields["VGENE_Reading_Frames: FR1,CDR1,FR2,CDR2,FR3,CDR3"]] = vGeneResults.readingFrameText["FR1"] + "," + vGeneResults.readingFrameText["CDR1"] + "," + vGeneResults.readingFrameText["FR2"] + "," + vGeneResults.readingFrameText["CDR2"] + "," + vGeneResults.readingFrameText["FR3"] + "," + vGeneResults.readingFrameText["CDR3"];
					//###summarize vgene resutls###/			
					vGermlineAlignment.AddUniqueLociToVectorString(uniqueLocusVector, currentLocusNum);
				}
				else{
					notes += "No alignment to VGene could be identified;";
				}
			}
			else{
				notes += "No alignment to VGene could be identified;";
			}
		}

		if (method["jgene"] != NONE && seqRead.seq.length()>minJAlignment){//Perform alignment to J GENES
			jStart = numVGermlineHits ==0 || vGeneResults.QE - 10 < 0 ? 0 : vGeneResults.QE - 10;
			dir = numVGermlineHits == 0 ? '\0' : dir;
			if (dir == '-'){
				locusSearch = vGermlineAlignment.NumberUniqueLoci() == 1 ? vGermlineAlignment.FirstLoci() : "";
				jGermlineAlignment.ResetQueryInfo(seqRead, 0, seqRead.seq.length() - jStart); //jstart reverse to start position along the reverse complement
				numJClusterHits = jGermlineAlignment.FindBestAlignmentDiagonal(-1,locusSearch);
			}
			else if (dir == '+'){
				locusSearch = vGermlineAlignment.NumberUniqueLoci() == 1 ? vGermlineAlignment.FirstLoci() : "";
				jGermlineAlignment.ResetQueryInfo(seqRead, jStart, seqRead.seq.length());
				numJClusterHits = jGermlineAlignment.FindBestAlignmentDiagonal(1,locusSearch);
			}
			else{
				jGermlineAlignment.ResetQueryInfo(seqRead, 0, seqRead.seq.length());
				numJClusterHits = jGermlineAlignment.FindBestAlignmentDiagonal();
			}
			if (numJClusterHits > 0){
				jGermlineAlignment.OptimizeJGeneClusterAlignmentWithGaps();
				numJGermlineHits = jGermlineAlignment.AlignClusterResultsToGermline();
				if (numJGermlineHits > 0){
					jGermlineAlignment.SummarizeJGeneResults(jGeneResults);


					if (numVGermlineHits == 0){
						abstart = jGeneResults.QS;//no v gene alignmetnw as found, so this is where it starts
						codonStart = jGeneResults.codonStart;
						startingAnnotation = "FR4";
					}
					endingAnnotation = "FR4";
					abend = jGeneResults.QE;
					dir = jGeneResults.direction;
					results_to_text[output_file_fields["Direction"]] = jGeneResults.direction;
					results_to_text[output_file_fields["Top_J-Gene_Hits"]] = jGeneResults.genes;
					results_to_text[output_file_fields["J-Gene_Alignment_Scores"]] = jGeneResults.scores;


					results_to_text[output_file_fields["JGENE: Query_Start"]] = std::to_string(jGeneResults.QS);
					results_to_text[output_file_fields["JGENE: Query_End"]] = std::to_string(jGeneResults.QE);
					results_to_text[output_file_fields["JGENE: Germline_Start"]] = std::to_string(jGeneResults.GS);
					results_to_text[output_file_fields["JGENE: Germline_End"]] = std::to_string(jGeneResults.GE);

					results_to_text[output_file_fields["JGENE: Alignment_Sequence_Query"]] = jGeneResults.queryAlSeq;
					results_to_text[output_file_fields["JGENE: Alignment_Sequence_Germline"]] = jGeneResults.germlineAlSeq;

					results_to_text[output_file_fields["JGENE: Total_Matches"]] = std::to_string(jGeneResults.totalMatch);
					results_to_text[output_file_fields["JGENE: Total_Mismatches"]] = std::to_string(jGeneResults.totalMismatch);
					results_to_text[output_file_fields["JGENE: Total_Indel"]] = std::to_string(jGeneResults.totalIndel);
					jGermlineAlignment.AddUniqueLociToVectorString(uniqueLocusVector, currentLocusNum);
				}
				else{
					notes += "No alignment to JGene could be identified;";
				}
			}
			else{
				notes += "No alignment to JGene could be identified;";
			}
		}

		results_to_text[output_file_fields["Header"]] = seqRead.seqHeader;
		results_to_text[output_file_fields["Sequence"]] = seqRead.seq;
		results_to_text[output_file_fields["Notes"]] = notes;
		results_to_text[output_file_fields["Sequence quality"]] = seqRead.quality;

		
		if (numJGermlineHits + numVGermlineHits > 0){
			if (dir == '-')
				seqRead.seq = dnafunctions::ReverseComplement(seqRead.seq);
			codonFrame = codonStart % 3 + 1;
			
			if (vGeneResults.cdr3start > -1){
				results_to_text[output_file_fields["CDR3:FR4_SEQUENCE.NT"]] = seqRead.seq.substr(vGeneResults.cdr3start, abend - vGeneResults.cdr3start + 1);
				results_to_text[output_file_fields["CDR3:FR4_SEQUENCE.AA"]] = dnafunctions::TranslateSeq(results_to_text[output_file_fields["CDR3:FR4_SEQUENCE.NT"]], vGeneResults.readingFrame["CDR3"], vGeneResults.cdr3start); // seqRead.seq.substr(abstart, abend - abstart + 1);
			}

			results_to_text[output_file_fields["Strand_Corrected_Sequence"]] = seqRead.seq;
			results_to_text[output_file_fields["Full_Length_Sequence.NT"]] = seqRead.seq.substr(abstart, abend - abstart + 1);						
			results_to_text[output_file_fields["Full_Length_Sequence: Start::End"]] = std::to_string(abstart) + "::" + std::to_string(abend);
			results_to_text[output_file_fields["Codon_Start"]] = std::to_string(codonStart);
			results_to_text[output_file_fields["Codon_Frame"]] = std::to_string(codonFrame);
			results_to_text[output_file_fields["5_Prime_Annotation"]] = startingAnnotation;
			results_to_text[output_file_fields["3_Prime_Annotation"]] = endingAnnotation;
			results_to_text[output_file_fields["Query_Start"]] = std::to_string(abstart);
			results_to_text[output_file_fields["Query_End"]] = std::to_string(abend);

			results_to_text[output_file_fields["Full_Length_Sequence.AA"]] = dnafunctions::TranslateSeq(results_to_text[output_file_fields["Full_Length_Sequence.NT"]], codonFrame,abstart); // seqRead.seq.substr(abstart, abend - abstart + 1);
			

			for (int k = 0; k < currentLocusNum - 1;k++)
				results_to_text[output_file_fields["Locus"]] += uniqueLocusVector[k] + ",";
			results_to_text[output_file_fields["Locus"]] += uniqueLocusVector[currentLocusNum - 1];			
		}
		WriteRow(outfile, results_to_text, headers.size());

		///++seqAnalyzed;
	}
	/*******************END OF ANALYSIS PER SEQUENCE *************/
	clock_t end_total_time = clock();
	printf("couts: %i", counts);
	printf("hits: %i", totalHits);

	//printf("First clock: %i", int(t));
	//printf("Second clock: %i", int(vGermlineAlignment.ReturnTime()));
	printf("FFT OUTSIDE FUNCTION: %i\n", int(t));
	printf("FFT INSIDE  FUNCTION: %i\n", int(vGermlineAlignment.ReturnTime()));
	printf("OPTIMIZE OUTSIDE  FUNCTION: %i\n", int(tna));
	printf("OPTIMIZE INSIDE  FUNCTION: %i\n", int(vGermlineAlignment.ReturnTimeOp()));
	printf("Total Time: %i\n", int(end_total_time - total_time));

	inFile.close();
	outfile.close();
}

vector<string> WriteHeaderRow(ofstream & outfile){
	vector<string> fields = { "Header", "Sequence","Sequence quality", "Direction", "Strand_Corrected_Sequence", "Query_Start", "Query_End", "Codon_Start", "Codon_Frame", "5_Prime_Annotation", "3_Prime_Annotation", "Full_Length_Sequence: Start::End", "Locus", "Top_V-Gene_Hits", "V-Gene_Alignment_Scores", "Top_J-Gene_Hits", "J-Gene_Alignment_Scores", "Full_Length_Sequence.NT", "Full_Length_Sequence.AA", "VGENE_Reading_Frames: FR1,CDR1,FR2,CDR2,FR3,CDR3", "FR1_Sequence.NT", "CDR1_Sequence.NT", "FR2_Sequence.NT", "CDR2_Sequence.NT", "FR3_Sequence.NT", "CDR3:FR4_SEQUENCE.NT", "FR1_Sequence.AA", "CDR1_Sequence.AA", "FR2_Sequence.AA", "CDR2_Sequence.AA", "FR3_Sequence.AA", "CDR3:FR4_SEQUENCE.AA", "VGENE: Total_Matches", "VGENE: Total_Mismatches", "VGENE: Total_Indel", "JGENE: Total_Matches", "JGENE: Total_Mismatches", "JGENE: Total_Indel", "VGENE_Matches: FR1,CDR1,FR2,CDR2,FR3,CDR3", "VGENE_Mismatches: FR1,CDR1,FR2,CDR2,FR3,CDR3", "VGENE_Indels: FR1,CDR1,FR2,CDR2,FR3,CDR3", "VGENE: Query_Start", "VGENE: Query_End", "VGENE: Query_FR1_Start::End", "VGENE: Query_CDR1_Start::End", "VGENE: Query_FR2_Start::End", "VGENE: Query_CDR2_Start::End", "VGENE: Query_FR3_Start::End", "VGENE: Germline_Start", "VGENE: Germline_End", "JGENE: Query_Start", "JGENE: Query_End", "JGENE: Germline_Start", "JGENE: Germline_End", "VGENE: Alignment_Sequence_Query", "VGENE: Alignment_Sequence_Germline", "VGENE: Alignment_FR1_Start::End", "VGENE: Alignment_CDR1_Start::End", "VGENE: Alignment_FR2_Start::End", "VGENE: Alignment_CDR2_Start::End", "VGENE: Alignment_FR3_Start::End", "JGENE: Alignment_Sequence_Query", "JGENE: Alignment_Sequence_Germline", "VGENE: Germline_Length", "JGENE: Germline_Length", "Notes" };

	for (int i = 0; i < fields.size() - 1; i++){
		outfile << fields[i] << "\t";
		output_file_fields[fields[i]] = i;
	}
 
	outfile << fields[fields.size() - 1] << "\n";
	output_file_fields[fields[fields.size() - 1]] = fields.size() - 1;

	return fields;
}

void WriteRow(ofstream & outfile, vector<string> fields, int numFields){
	for (int i = 0; i < numFields - 1; i++){
		outfile << fields[i] << "\t";
	}
	outfile << fields[numFields - 1] << "\n";
}

void SummarizeData(vector<string> & fields, const vector<structvars::Germline> & germline_set, const vector<vector<double>> & germline_scores, const int & numGenes){
	fields[3] = "";
	fields[4] = "";

	for (int i = 0; i < numGenes; i++){
		fields[3] += germline_set[germline_scores[i][1]].genename + ",";
		fields[4] += to_string(germline_scores[i][0]) + ",";
	}
}

//This function will evaluate the parameters inputed from the user
void EvaluateParameters(int argc, char *argv[]){
	char var;
	int currentVar = 1;
	string parameter, headVal;
	//	int parameterVarInt;
	FILE *pFile;
	string filename; //for reading filenames/strings
	std::istringstream ss;
	string token; //use this for istringsstream, see below
	int numDatabase;
	//printf("%i \n",argc);
	//for (int k =0;k<argc;k++)
	//	printf("%s \n",argv[k]);
	if (argc < 4){//check input arguments
		if (string(argv[1]) == "--defaults"){
			PrintDefaultParameters();
			exit(1);
		}
		else if (string(argv[1]) == "--version")
		{
			printf("Version 0.6");
			exit(1);
		}
		else{

			printf("Error: You must enter the location of the sequence file, and the location of at least one database to compare to\n");
			exit(EXIT_FAILURE);
		}
	}
	else if (argc % 2 != 0){//check if even number of arguments
		printf("Error entering parameters\n");
		exit(EXIT_FAILURE);
	}
	pFile = fopen(argv[1], "r");
	if (!pFile){//check if input file name exists
		printf("Error opening file\n");
		exit(EXIT_FAILURE);
	}
	fclose(pFile);
	std::printf("Evaluating Parameters\n");
	variables_used.fileinputName = argv[1];//input file name
	variables_used.fileoutputName = variables_used.fileinputName + ".fft.annotation";
	while (currentVar < argc - 1){
		currentVar++;
		parameter = string(argv[currentVar]);
		var = ::tolower(parameter[1]);// ::tolower(argv[currentVar][1]);
		if (argv[currentVar][0] != '-'){
			printf("Incorrect parameter\n");
			exit(EXIT_FAILURE);
		}
		cout << argv[currentVar] << endl;
		switch (var)
		{

		case 'v': //location of the vgermline database
			currentVar++;
			if (strlen(argv[currentVar - 1]) == 2){
				filename = string(argv[currentVar]);
				ss.clear();
				ss.str(filename);
				numDatabase = 0;
				//input each of the filenames delimited by a ','
				while (getline(ss, token, ',')){
					try{
						if (token != ""){
							readfiles::RemoveTrailingSpaces(token); //remove all trailing white spaces
							//make sure teh files they passed in exist
							pFile = fopen(token.c_str(), "r");
							if (!pFile){//check if input file name exists
								printf("V Germline database file location does not exist: %s", token.c_str());
								exit(EXIT_FAILURE);
							}
							fclose(pFile);
							variables_used.query_algn_settings['V'].fileGermline[numDatabase] = token; //the first index MUST be the location of the list of germline sequences.  the second index MUST be the location of the clustered germline seuqences.
							numDatabase++;
						}
					}
					catch (std::exception&e){
						printf("The V Gene Database must be defined by two file locations: the first location is the list of germline sequences, and the second file refers to the clustered germline sequences");
						exit(EXIT_FAILURE);
					}
				}
				//error out if we did not read 2 files (index starts at 0)
				if (numDatabase < 1 || variables_used.query_algn_settings['V'].fileGermline[0] == "")
				{
					printf("The V Gene Database file location must be defined");
					exit(EXIT_FAILURE);
				}
				//define whether user desires vgene data
				variables_used.containsVGermline = true;
			}
			else{
				EvaluateAlgnParameters('V', parameter.substr(2), string(argv[currentVar]));
			}
			break;
		case 'j': //location of the jgermline database
			currentVar++;
			//define whether the user desires j gene data
			if (strlen(argv[currentVar - 1]) == 2){
				variables_used.containsJGermline = true;
				filename = string(argv[currentVar]);
				ss.clear();
				ss.str(filename);
				numDatabase = 0;
				//input each of the filenames delimited by a ','
				while (getline(ss, token, ',')){
					try{
						//token.erase(remove_if(token.begin(), token.end(), isspace), token.end());
						readfiles::RemoveTrailingSpaces(token); //remove trailling white spaces
						//make sure teh files they passed in exist
						pFile = fopen(token.c_str(), "r");
						if (!pFile){//check if input file name exists
							printf("J Germline database file location does not exist: %s", token.c_str());
							exit(EXIT_FAILURE);
						}
						fclose(pFile);
						variables_used.query_algn_settings['J'].fileGermline[numDatabase] = token; //the first index MUST be the location of the list of germline sequences.  the second index MUST be the location of the clustered germline seuqences.
						numDatabase++;
					}
					catch (std::exception&e){
						printf("The J Gene Database must be defined by two file locations: the first location is the list of germline sequences, and the second file refers to the clustered germline sequences");
						exit(EXIT_FAILURE);
					}
				}
				//error out if we did not read 2 files (index starts at 0)
				if (numDatabase < 1 || variables_used.query_algn_settings['J'].fileGermline[0] == "")
				{
					printf("The J Gene Database file location must be provided");
					exit(EXIT_FAILURE);
				}
			}
			else{
				EvaluateAlgnParameters('J', parameter.substr(2), string(argv[currentVar]));
			}
			break;
		case 'd': //location of the jgermline database
			currentVar++;
			//define whether the user desires j gene data
			if (strlen(argv[currentVar - 1]) == 2){
				variables_used.containsJGermline = true;

				filename = string(argv[currentVar]);
				ss.clear();
				ss.str(filename);
				numDatabase = 0;
				//input each of the filenames delimited by a ','
				while (getline(ss, token, ',')){
					try{
						//token.erase(remove_if(token.begin(), token.end(), isspace), token.end());
						readfiles::RemoveTrailingSpaces(token); //remove trailling white spaces
						//make sure teh files they passed in exist
						pFile = fopen(token.c_str(), "r");
						if (!pFile){//check if input file name exists
							printf("J Germline database file location does not exist: %s", token.c_str());
							exit(EXIT_FAILURE);
						}
						fclose(pFile);
						variables_used.query_algn_settings['D'].fileGermline[numDatabase] = token; //the first index MUST be the location of the list of germline sequences.  the second index MUST be the location of the clustered germline seuqences.
						numDatabase++;
					}
					catch (std::exception&e){
						printf("The J Gene Database must be defined by two file locations: the first location is the list of germline sequences, and the second file refers to the clustered germline sequences");
						exit(EXIT_FAILURE);
					}
				}
				//error out if we did not read 2 files (index starts at 0)
				if (numDatabase < 1 || variables_used.query_algn_settings['D'].fileGermline[0] == "")
				{
					printf("The J Gene Database file location must be provided");
					exit(EXIT_FAILURE);
				}
			}
			else{
				EvaluateAlgnParameters('D', parameter.substr(2), string(argv[currentVar]));
			}
			break;
		case 'o': //location of the file output name
			currentVar++;
			variables_used.fileoutputName = string(argv[currentVar]);
			break;
		case 'i': //input file format
			currentVar++;
			try{
				variables_used.inputFileFormat = string(argv[currentVar]);
				for (int k = 0; k < variables_used.inputFileFormat.length(); k++){
					variables_used.inputFileFormat[k] = ::toupper(variables_used.inputFileFormat[k]);
				}

				if (variables_used.inputFileFormat != "FASTA" && variables_used.inputFileFormat != "TAB" && variables_used.inputFileFormat != "FASTQ"){
					printf("The only allowed values for file format are 'TAB' for text tab delimited file and 'FASTA' for a fasta file");
					exit(EXIT_FAILURE);
				};
			}
			catch (std::exception& e)
			{
				printf("The only allowed values for file format are 'TAB' for text tab delimited file and 'FASTA' for a fasta file");
				exit(EXIT_FAILURE);
			}
			break;
		case 'h': // pass in variable for whether or not to ignore the first line in a tab file
			currentVar++;
			headVal = string(argv[currentVar]);
			for (int k = 0; k < headVal.length(); k++){
				headVal[k] = ::tolower(headVal[k]);
			}
			if (headVal != "false" && headVal != "true")
			{
				printf("The only allowed values for ignore header parameter (-h) are 'true' or 'false'");
				exit(EXIT_FAILURE);
			}
			variables_used.ignoreHeader = headVal == "false" ? false : true;
			break;
		default:
			currentVar++;
			EvaluateAlgnParameters('\0', parameter.substr(1), string(argv[currentVar]));
			break;
		}
	}
	//check if either of the databases were included. at least one must be defined.
	//that is are they both false....
	if (!(variables_used.containsJGermline || variables_used.containsVGermline)){
		printf("Error: You must enter the location of at least one germline database\n");
		exit(EXIT_FAILURE);
	}
}

void PrintDefaultParameters(bool to_file){
	string output_text = "";
	structvars::AlignmentProgramSettings defaultalgnsettings;

	output_text += "-v\tnone\tLocation of V Germline Database\n-d\tnone\tLocation of D Germline Database\n-j\tnone\tLocation of D Germline Database\n";
	output_text += "-h\ttrue\tQuery file includes header row\n-i\tFASTA\tInput file type\n-o\tnone\toutput filename\n";

	typedef std::map<std::string, settings_description >::iterator it_type;
	for (it_type iterator = interpreter.begin(); iterator != interpreter.end(); iterator++) {
		UserSettings eval = iterator->second.parameter;
		string value;

		switch (eval){
		case allowed_gaps:
			value = std::to_string(defaultalgnsettings.fftParams.maxGap);
			break;
		case match_sw:
			value = std::to_string(defaultalgnsettings.swParams.matchScore);
			break;
		case mismatch_sw:
			value = std::to_string(defaultalgnsettings.swParams.mismatchScore);
			break;
		case gap_open_sw:
			value = std::to_string(defaultalgnsettings.swParams.swGapOpen);
			break;
		case gap_extend_sw:
			value = std::to_string(defaultalgnsettings.swParams.swExtendGap);
			break;
		case gap_extend_fft:
			value = std::to_string(defaultalgnsettings.fftParams.fftExtendGap);
			break;
		case gap_open_fft:
			value = std::to_string(defaultalgnsettings.fftParams.fftGapOpen);
			break;
		case sensitivity_fft:
			value = std::to_string(defaultalgnsettings.fftParams.sensitivity);
			break;
		case sensitivity_peptide:
			value = std::to_string(defaultalgnsettings.peptideParams.peptide_gap_sensitivity);
			break;
		case similar_clusters:
			value = std::to_string(defaultalgnsettings.fftParams.cluster_threshold);
			break;
		case peptide_len:
			value = std::to_string(defaultalgnsettings.peptideParams.peptide_len);
			break;
		case min_per_id_threshold:
			value = std::to_string(defaultalgnsettings.fftParams.scoreCutoff);
			break;
		case max_germline_hits:
			value = std::to_string(defaultalgnsettings.maxHits);
			break;
		case score_ratio:
			value = std::to_string(defaultalgnsettings.scoreFoldRatio);
			break;
		case num_above_score_ratio:
			value = std::to_string(defaultalgnsettings.numObsAboveFoldRatio);
			break;
		case cluster_germlines_cutoff:
			value = std::to_string(defaultalgnsettings.clusterGermlineCutoff);
			break;
		case group_clusters:
			value = defaultalgnsettings.group_into_clusters ? "true" : "false" ;
			break;
		};

		output_text += "-" + iterator->first + "\t" + value + "\t" + iterator->second.description + "\n";
		// iterator->first = key
		// iterator->second = value
		// Repeat if you also want to iterate through the second map.
	}

	if (to_file){
		ofstream outfile("defaultsettings_fftprogram.txt");
		outfile << output_text;
		outfile.close();
	}
	else{
		printf("\%s", output_text.c_str());
	}

}

void EvaluateAlgnParameters(char germline, string parameter, string value){
	UserSettings eval = interpreter[parameter].parameter;

	//Just a copy of the the variable interpreter above for reference
	/*
	interpreter["gap"] = allowed_gaps;
	interpreter["match_sw"] = match_sw;
	interpreter["mismatch_sw"] = mismatch_sw;
	interpreter["gap_open_sw"] = gap_open_sw;
	interpreter["gap_extend_sw"] = gap_extend_sw;
	interpreter["gap_open_fft"] = gap_open_fft;
	interpreter["gap_extend_fft"] = gap_extend_fft;
	interpreter["s_fft"] = sensitivity_fft;
	interpreter["s_pep"] = sensitivity_peptide;
	interpreter["similar_clusters"] = similar_clusters;
	interpreter["pep_len"] = peptide_len;
	interpreter["min_per_id"] = min_per_id_threshold;
	interpreter["num_hits"] = max_germline_hits;
	interpreter["ratio_cutoff"] = score_ratio;
	interpreter["times_above_ratio"] = num_above_score_ratio;
	interpreter["cluster_per_id"] = cluster_germlines_cutoff;
	interpreter["group_clusters"] = group_clusters;
	*//////////////////////////////////////////////////////////

	int temp_i;
	double temp_d;

	switch (eval){
	case allowed_gaps:
		temp_i = abs(GetIntegerValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].fftParams.maxGap = temp_i;
			variables_used.query_algn_settings['V'].swParams.maxGap = temp_i;
			variables_used.query_algn_settings['D'].fftParams.maxGap = temp_i;
			variables_used.query_algn_settings['D'].swParams.maxGap = temp_i;
			variables_used.query_algn_settings['J'].fftParams.maxGap = temp_i;
			variables_used.query_algn_settings['J'].swParams.maxGap = temp_i;
		}
		else{
			variables_used.query_algn_settings[germline].fftParams.maxGap = temp_i;
			variables_used.query_algn_settings[germline].swParams.maxGap = temp_i;
		}
		break;
	case match_sw:
		temp_i = abs(GetIntegerValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].swParams.matchScore = temp_i;
			variables_used.query_algn_settings['D'].swParams.matchScore = temp_i;
			variables_used.query_algn_settings['J'].swParams.matchScore = temp_i;
		}
		else{
			variables_used.query_algn_settings[germline].swParams.matchScore = temp_i;
		}
		break;
	case mismatch_sw:
		temp_i = abs(GetIntegerValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].swParams.mismatchScore = temp_i;
			variables_used.query_algn_settings['D'].swParams.mismatchScore = temp_i;
			variables_used.query_algn_settings['J'].swParams.mismatchScore = temp_i;
		}
		else{
			variables_used.query_algn_settings[germline].swParams.mismatchScore = temp_i;
		}
		break;
	case gap_open_sw:
		temp_d = abs(GetDecimalValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].swParams.swGapOpen = temp_d;
			variables_used.query_algn_settings['D'].swParams.swGapOpen = temp_d;
			variables_used.query_algn_settings['J'].swParams.swGapOpen = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].swParams.swGapOpen = temp_d;
		}
		break;
	case gap_extend_sw:
		temp_d = abs(GetDecimalValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].swParams.swExtendGap = temp_d;
			variables_used.query_algn_settings['D'].swParams.swExtendGap = temp_d;
			variables_used.query_algn_settings['J'].swParams.swExtendGap = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].swParams.swExtendGap = temp_d;
		}
		break;
	case gap_extend_fft:
		temp_d = abs(GetDecimalValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].fftParams.fftExtendGap = temp_d;
			variables_used.query_algn_settings['D'].fftParams.fftExtendGap = temp_d;
			variables_used.query_algn_settings['J'].fftParams.fftExtendGap = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].fftParams.fftExtendGap = temp_d;
		}
		break;
	case gap_open_fft:
		temp_d = abs(GetDecimalValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].fftParams.fftGapOpen = temp_d;
			variables_used.query_algn_settings['D'].fftParams.fftGapOpen = temp_d;
			variables_used.query_algn_settings['J'].fftParams.fftGapOpen = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].fftParams.fftGapOpen = temp_d;
		}
		break;
	case sensitivity_fft:
		temp_d = abs(GetDecimalValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].fftParams.sensitivity = temp_d;
			variables_used.query_algn_settings['D'].fftParams.sensitivity = temp_d;
			variables_used.query_algn_settings['J'].fftParams.sensitivity = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].fftParams.sensitivity = temp_d;
		}
		break;
	case sensitivity_peptide:
		temp_d = abs(GetDecimalValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].peptideParams.peptide_gap_sensitivity = temp_d;
			variables_used.query_algn_settings['D'].peptideParams.peptide_gap_sensitivity = temp_d;
			variables_used.query_algn_settings['J'].peptideParams.peptide_gap_sensitivity = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].peptideParams.peptide_gap_sensitivity = temp_d;
		}
		break;
	case similar_clusters:
		temp_d = abs(GetDecimalValue(value, true));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].fftParams.cluster_threshold = temp_d;
			variables_used.query_algn_settings['D'].fftParams.cluster_threshold = temp_d;
			variables_used.query_algn_settings['J'].fftParams.cluster_threshold = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].fftParams.cluster_threshold = temp_d;
		}
		break;
	case peptide_len:
		temp_i = abs(GetIntegerValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].peptideParams.peptide_len = temp_i;
			variables_used.query_algn_settings['D'].peptideParams.peptide_len = temp_i;
			variables_used.query_algn_settings['J'].peptideParams.peptide_len = temp_i;
		}
		else{
			variables_used.query_algn_settings[germline].peptideParams.peptide_len = temp_i;
		}
		break;
	case min_per_id_threshold:
		temp_d = abs(GetDecimalValue(value, true));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].fftParams.scoreCutoff = temp_d;
			variables_used.query_algn_settings['D'].fftParams.scoreCutoff = temp_d;
			variables_used.query_algn_settings['J'].fftParams.scoreCutoff = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].fftParams.scoreCutoff = temp_d;
		}
		break;
	case max_germline_hits:
		temp_i = abs(GetIntegerValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].maxHits = temp_i;
			variables_used.query_algn_settings['D'].maxHits = temp_i;
			variables_used.query_algn_settings['J'].maxHits = temp_i;
		}
		else{
			variables_used.query_algn_settings[germline].maxHits = temp_i;
		}
		break;
	case score_ratio:
		temp_d = abs(GetDecimalValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].scoreFoldRatio = temp_d;
			variables_used.query_algn_settings['D'].scoreFoldRatio = temp_d;
			variables_used.query_algn_settings['J'].scoreFoldRatio = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].scoreFoldRatio = temp_d;
		}
		break;
	case num_above_score_ratio:
		temp_i = abs(GetIntegerValue(value, false));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].numObsAboveFoldRatio = temp_i;
			variables_used.query_algn_settings['D'].numObsAboveFoldRatio = temp_i;
			variables_used.query_algn_settings['J'].numObsAboveFoldRatio = temp_i;
		}
		else{
			variables_used.query_algn_settings[germline].numObsAboveFoldRatio = temp_i;
		}
		break;
	case cluster_germlines_cutoff:
		temp_d = abs(GetDecimalValue(value, true));
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].clusterGermlineCutoff = temp_d;
			variables_used.query_algn_settings['D'].clusterGermlineCutoff = temp_d;
			variables_used.query_algn_settings['J'].clusterGermlineCutoff = temp_d;
		}
		else{
			variables_used.query_algn_settings[germline].clusterGermlineCutoff = temp_d;
		}
		break;
	case group_clusters:
		for (int k = 0; k < value.length(); k++){
			value[k] = ::tolower(value[k]);
		}
		if (value != "false" && value != "true")
		{
			printf("The only allowed values for %s are 'true' or 'false'", parameter.c_str());
			exit(EXIT_FAILURE);
		}
		if (germline == '\0'){
			variables_used.query_algn_settings['V'].group_into_clusters = value == "false" ? false : true;
			variables_used.query_algn_settings['D'].group_into_clusters = value == "false" ? false : true;
			variables_used.query_algn_settings['J'].group_into_clusters = value == "false" ? false : true;
		}
		else{
			variables_used.query_algn_settings[germline].group_into_clusters = value == "false" ? false : true;
		}
		break;
	}

}

double GetDecimalValue(string value, bool between_0_1){
	double temp_d;
	string added_error_string = between_0_1 ? "from 0 to 1" : "";

	try{
		temp_d = std::stod(value);
	}
	catch (std::exception& e)
	{
		printf("%s must be a decimal number %s", value.c_str(), added_error_string.c_str());
		exit(EXIT_FAILURE);
	}

	if (between_0_1){
		if (temp_d < 0 || temp_d>1){
			printf("%s must be a decimal number %s", value.c_str(), added_error_string.c_str());
			exit(EXIT_FAILURE);
		}
	}
	return temp_d;
}
int GetIntegerValue(string value, bool between_0_1){
	int temp_i;
	string added_error_string = between_0_1 ? "from 0 to 1" : "";

	try{
		temp_i = std::stoi(value);
	}
	catch (std::exception& e)
	{
		printf("%s must be an integer %s", value.c_str(), added_error_string.c_str());
		exit(EXIT_FAILURE);
	}
	if (between_0_1){		
		if (temp_i < 0 || temp_i>1){			
			printf("%s must be an integer number %s", value.c_str(), added_error_string.c_str());		
			exit(EXIT_FAILURE);		
		}	
	}	
	return temp_i;
}
