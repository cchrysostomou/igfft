#define _CRT_SECURE_NO_DEPRECATE
#include "ReadFiles.h"
using namespace std;
using namespace readfiles;
using namespace structvars;
FASTQinfo readfiles::readIndividualSequence(std::ifstream & inFile)
{
	string header;
	string sequence;
	int stLen, g, asciVal;
	FASTQinfo fastqSequence;
	//READ IN HEADER SEQUENCE
	getline(inFile, header); //read the header until end of line
	if (header[0] == '>')//if we took in the '>' symbol, remove it
		header = header.substr(1, header.length() - 1);
	//READ IN ACTUAL SEQUENCE
	getline(inFile, sequence, '>');
	//now we ahve to clean up "sequence" and remove white spaces and "new line" characters
	stLen = sequence.length();
	g = 0;
	while (g < stLen){
		asciVal = (int)sequence[g];
		if (!((asciVal >= 97 && asciVal <= 122) || (asciVal >= 65 && asciVal <= 90))){ //if it is NOT a letter
			sequence.erase(g, 1); //erase the character
			stLen = sequence.length(); //get the new character length
			g--; //move the position back one because we have just erased a letter
		}
		++g;
	}
	RemoveTrailingSpaces(sequence);
	RemoveTrailingSpaces(header);
	std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper); //convert string to uppercase
	fastqSequence.seq = sequence;
	fastqSequence.seqHeader = header;
	return fastqSequence;
}

FASTAinfo readfiles::readIndividualSequenceToFASTA(std::ifstream & inFile)
{
	string header;
	string sequence;
	int stLen, g, asciVal;
	FASTAinfo fastaSequence;
	//READ IN HEADER SEQUENCE
	getline(inFile, header); //read the header until end of line
	if (header[0] == '>')//if we took in the '>' symbol, remove it
		header = header.substr(1, header.length() - 1);
	//READ IN ACTUAL SEQUENCE
	getline(inFile, sequence, '>');
	//now we ahve to clean up "sequence" and remove white spaces and "new line" characters
	stLen = sequence.length();
	g = 0;
	while (g < stLen){
		asciVal = (int)sequence[g];
		if (!((asciVal >= 97 && asciVal <= 122) || (asciVal >= 65 && asciVal <= 90))){ //if it is NOT a letter
			sequence.erase(g, 1); //erase the character
			stLen = sequence.length(); //get the new character length
			g--; //move the position back one because we have just erased a letter
		}
		++g;
	}
	RemoveTrailingSpaces(sequence);
	RemoveTrailingSpaces(header);
	std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper); //convert string to uppercase
	fastaSequence.seq = sequence;
	fastaSequence.seqHeader = header;
	return fastaSequence;
}

FASTQinfo readfiles::readIndividualFASTQSequence(std::ifstream & inFile)
{
	//structure of FASTQ: 
		//@header
		//sequence 
		//+
		//quality

	string header;
	string sequence;
	string quality;
	string temp;
	int stLen, g, asciVal,qualLen=0;
	FASTQinfo fastqSequence;
	//READ IN HEADER SEQUENCE
	getline(inFile, header); //read the header until end of line
	
	if (header.length() == 0){//no line found...skip result
		return fastqSequence;
	}
	else if (header[0] == '@')//if we took in the '@' symbol, remove it
		header = header.substr(1, header.length() - 1);
	else{
		cout << "ERROR QUALITY LINE DOES NOT START WITH '@'";
		return fastqSequence;
	}
	//READ IN ACTUAL SEQUENCE //Go up until the '+' sign 
	getline(inFile, sequence, '+');
	//now we ahve to clean up "sequence" and remove white spaces and "new line" characters
	stLen = sequence.length();
	g = 0;
	while (g < stLen){
		asciVal = (int)sequence[g];
		if (!((asciVal >= 97 && asciVal <= 122) || (asciVal >= 65 && asciVal <= 90))){ //if it is NOT a letter
			sequence.erase(g, 1); //erase the character
			stLen = sequence.length(); //get the new character length
			g--; //move the position back one because we have just erased a letter
		}
		++g;
	}
	RemoveTrailingSpaces(sequence);

	stLen = sequence.length(); //this is the length of the sequence. Use this length to read quality score
	while (qualLen < stLen){
		//read in quality value
		getline(inFile, temp);
		quality.append(temp);
		RemoveTrailingSpaces(quality);
		qualLen = quality.length();
	}

	RemoveTrailingSpaces(header);
	if (qualLen != stLen)
		cout << "ERROR: LENGTH OF QUALITY STRING DOES NOT MATCH SEQUENCE LENGTH" << endl;

	std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper); //convert string to uppercase
	fastqSequence.seq = sequence;
	fastqSequence.seqHeader = header;
	fastqSequence.quality = quality;
	return fastqSequence;
}

//Same as aove function, excpet for a tab delimited sequence
FASTQinfo readfiles::readIndividualTabSequence(std::ifstream & inFile){
	string header, line,quality;
	string sequence;
	FASTQinfo fastqSequence;
	std::istringstream ss;
	//READ IN HEADER SEQUENCE
	getline(inFile, line); //read the header until end of line
		
	ss.str(line);
	if(line==""){
		fastqSequence.seqHeader="Noneerror";
		fastqSequence.seq="";
		return fastqSequence;
	}
	
	getline(ss, header, '\t');
	getline(ss, sequence, '\t');
	getline(ss, quality, '\t');
	RemoveTrailingSpaces(header);
	RemoveTrailingSpaces(sequence);	
	RemoveTrailingSpaces(quality);
	std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper); //convert string to uppercase
	fastqSequence.seq = sequence;	
	fastqSequence.quality = quality;
	fastqSequence.seqHeader = header;
	return fastqSequence;
}

//function read an entire fasta file and add it to a vector
Fileinfo readfiles::readFASTAToVector(const char* filename, vector<FASTQinfo> &sequence)
{
	std::ifstream inFile;
	inFile.open(filename); //open the FASTA file
	FASTQinfo fastaSequence;
	Fileinfo results;
	string temp;
	while (true){
		getline(inFile, temp); //read each line 		
		if (temp[0] == '>'){//START OF FASTA FILE!!
			inFile.putback('\n'); //PUT BACK EOL CHARACTER!!!
			for (int i = temp.length()-1; i >= 0;i--)		
				inFile.putback(temp[i]);
			break;
		}
	}
	
	
	int numSeq = 0, minLen = 1000, maxLen = 0, totalLen=0;
	while (!inFile.eof()) {
		fastaSequence = readIndividualSequence(inFile);
		if (fastaSequence.seq != "" && fastaSequence.seqHeader != ""){ //just make sure we actually read a proper line of sequence and sequence header
			sequence.push_back(fastaSequence);
			if (fastaSequence.seq.length() < minLen)
				minLen = fastaSequence.seq.length();
			if (fastaSequence.seq.length()>maxLen)
				maxLen = fastaSequence.seq.length();
			AddSeqLenToVecDistribution(fastaSequence.seq.length(), results.lengthDistribution);
			totalLen += fastaSequence.seq.length();
			++numSeq; //add one to the number of new sequences
		}
	}
	inFile.close();
	results.avgSeqLen = (totalLen / numSeq) + 1;
	results.maxSeqLen = maxLen;
	results.minSeqLen = minLen;
	results.numSeqs = numSeq;
	printf("A total of %i sequences are present in the query file\n", numSeq);
	sort(results.lengthDistribution.begin(), results.lengthDistribution.end(), [](const vector<int> & a, const vector<int> & b){ return (a[3] > b[3]); });
	return results;
}

//function read an entire fasta file and add it to a vector
Fileinfo readfiles::readFASTQToVector(const char* filename, vector<FASTQinfo> &sequence,char filetype)
{
	std::ifstream inFile;
	string temp;
	inFile.open(filename); //open the FASTQ file
	FASTQinfo fastqSequence;
	Fileinfo results;
	int numSeq = 0, minLen = 1000, maxLen = 0, totalLen = 0;
	
	while (true){
		getline(inFile, temp); //read each line 		
		if ((filetype == 'A' && temp[0] == '>') || (filetype == 'Q' && temp[0] == '@')){ //START of FASTA file  or FAASTQ file 
			inFile.putback('\n');//PUT BACK EOL character
			for (int i = temp.length()-1; i >= 0; i--)
				inFile.putback(temp[i]);
			break;
		}		
	}	

	while (!inFile.eof()) {
		if (filetype == 'A')
			fastqSequence = readIndividualSequence(inFile);
		else
			fastqSequence = readIndividualFASTQSequence(inFile);
		if (fastqSequence.seq != "" && fastqSequence.seqHeader != ""){ //just make sure we actually read a proper line of sequence and sequence header
			sequence.push_back(fastqSequence);
			if (fastqSequence.seq.length() < minLen)
				minLen = fastqSequence.seq.length();
			if (fastqSequence.seq.length()>maxLen)
				maxLen = fastqSequence.seq.length();
			AddSeqLenToVecDistribution(fastqSequence.seq.length(), results.lengthDistribution);
			totalLen += fastqSequence.seq.length();
			++numSeq; //add one to the number of new sequences
		}
	}
	inFile.close();
	results.avgSeqLen = (totalLen / numSeq) + 1;
	results.maxSeqLen = maxLen;
	results.minSeqLen = minLen;
	results.numSeqs = numSeq;
	printf("A total of %i sequences are present in the query file\n", numSeq);
	sort(results.lengthDistribution.begin(), results.lengthDistribution.end(), [](const vector<int> & a, const vector<int> & b){ return (a[3] > b[3]); });
	return results;
}

//function read an entire tab file and adds it to a vector
Fileinfo readfiles::readTABToVector(const char* filename, vector<FASTQinfo> &sequence)
{
	std::ifstream inFile;
	inFile.open(filename); //open the FASTA file
	FASTQinfo fastqSequence;
	Fileinfo results;
	int numSeq = 0, minLen = 1000, maxLen = 0,totalLen = 0;
	while (!inFile.eof()) {
		fastqSequence = readIndividualTabSequence(inFile);
		if (fastqSequence.seq != "" && fastqSequence.seqHeader != ""){//just make sure we actually read a proper line of sequence and sequence header
			sequence.push_back(fastqSequence);
			if (fastqSequence.seq.length() < minLen)
				minLen = fastqSequence.seq.length();
			if (fastqSequence.seq.length()>maxLen)
				maxLen = fastqSequence.seq.length();
			AddSeqLenToVecDistribution(fastqSequence.seq.length(), results.lengthDistribution);
			totalLen += fastqSequence.seq.length();
			++numSeq; //add one to the number of new sequences
		}
	}
	inFile.close();
	results.avgSeqLen = (totalLen / numSeq) + 1;
	results.maxSeqLen = maxLen;
	results.minSeqLen = minLen;
	results.numSeqs = numSeq;
	printf("A total of %i sequences are present in the query file\n", numSeq);
	sort(results.lengthDistribution.begin(), results.lengthDistribution.end(), [](const vector<int> & a, const vector<int> & b){ return (a[3] > b[3]); });
	return results;
}

//convert a fasta file to a tab delimited file
structvars::Fileinfo readfiles::convertFASTAtoTAB(string & inputFileName){
	structvars::Fileinfo data;
	data.numSeqs = 0;
	string appendData = ".convert_tab.txt";
	string outputfilename = string(inputFileName);
	outputfilename.append(appendData);
	int minLen = 1000, maxLen = 0, totalLen = 0;
	ofstream outFile(outputfilename, ofstream::out);
	ifstream inFile;
	inFile.open(inputFileName.c_str());
	structvars::FASTAinfo tempFastaSequence;
	//skip all lines at beginnning of FASTA file (only start where file contains '<')
	string temp;
	while (true){
		getline(inFile, temp); //read each line 		
		if (temp[0] == '>'){//START OF FASTA FILE!!
			inFile.putback('\n');//PUT BACK EOL character		
			for (int i = temp.length()-1; i >= 0; i--)
				inFile.putback(temp[i]);
			break;
		}
	}
	while (!inFile.eof()){
		tempFastaSequence = readIndividualSequenceToFASTA(inFile);
		//if (tempFastaSequence.seqHeader != ""){//just make sure we actually read a proper line with a sequence header		
		if (tempFastaSequence.seq.length() > 0){
			WriteFASTASeqToTab(outFile, tempFastaSequence);
			if (tempFastaSequence.seq.length() < minLen && tempFastaSequence.seq.length() > 0)
				minLen = tempFastaSequence.seq.length();
			if (tempFastaSequence.seq.length() > maxLen)
				maxLen = tempFastaSequence.seq.length();
			AddSeqLenToVecDistribution(tempFastaSequence.seq.length(), data.lengthDistribution);
			totalLen += tempFastaSequence.seq.length();
			data.numSeqs += 1;
		}
		//}
		
	}
	inFile.close();
	data.avgSeqLen = (totalLen / data.numSeqs) + 1;
	data.minSeqLen = minLen;
	data.maxSeqLen = maxLen;
	printf("A total of %i sequences were found in the query file\n", data.numSeqs);
	sort(data.lengthDistribution.begin(), data.lengthDistribution.end(), [](const vector<int> & a, const vector<int> & b){ return (a[3] > b[3]); });
	inputFileName = outputfilename;
	return data;
}

//convert a fasta file to a tab delimited file
structvars::Fileinfo readfiles::convertFASTQtoTAB(string & inputFileName){
	structvars::Fileinfo data;
	data.numSeqs = 0;
	string appendData = ".convert_tab.txt";
	string outputfilename = string(inputFileName);
	outputfilename.append(appendData);
	int minLen = 1000, maxLen = 0, totalLen = 0;
	ofstream outFile(outputfilename, ofstream::out);
	ifstream inFile;
	inFile.open(inputFileName.c_str());
	structvars::FASTQinfo tempFastqSequence;
	//skip all lines at beginnning of FASTA file (only start where file contains '<')
	string temp;
	while (true){
		getline(inFile, temp); //read each line 		
		if (temp[0] == '@'){//START OF FASTQ FILE!!
			inFile.putback('\n');//PUT BACK EOL CHARACTER! 
			for (int i = temp.length()-1; i >= 0; i--)
				inFile.putback(temp[i]);
			break;
		}
	}
	
	while (!inFile.eof()){
		
		tempFastqSequence = readIndividualFASTQSequence(inFile);		
		if (tempFastqSequence.seq.length() > 0){
			//if (tempFastaSequence.seqHeader != ""){//just make sure we actually read a proper line with a sequence header
			WriteFASTQSeqToTab(outFile, tempFastqSequence);			
			if (tempFastqSequence.seq.length() < minLen && tempFastqSequence.seq.length() > 0)
				minLen = tempFastqSequence.seq.length();
			if (tempFastqSequence.seq.length() > maxLen)
				maxLen = tempFastqSequence.seq.length();
			AddSeqLenToVecDistribution(tempFastqSequence.seq.length(), data.lengthDistribution);
			totalLen += tempFastqSequence.seq.length();
			data.numSeqs += 1;
		}		
		
	}

	inFile.close();
	data.avgSeqLen = (totalLen / data.numSeqs) + 1;
	data.minSeqLen = minLen;
	data.maxSeqLen = maxLen;
	printf("A total of %i sequences were found in the query file\n", data.numSeqs);
	sort(data.lengthDistribution.begin(), data.lengthDistribution.end(), [](const vector<int> & a, const vector<int> & b){ return (a[3] > b[3]); });
	inputFileName = outputfilename;
	return data;
}

void readfiles::WriteFASTASeqToTab(ofstream & outFile, FASTQinfo fastaseq){
	string line;
	line = fastaseq.seqHeader + "\t" + fastaseq.seq + "\t\n";
	outFile << line;
}

void readfiles::WriteFASTQSeqToTab(ofstream & outFile, FASTQinfo fastqseq){
	string line;
	line = fastqseq.seqHeader + "\t" + fastqseq.seq + "\t"+fastqseq.quality+"\n";
	outFile << line;
}

void readfiles::WriteFASTASeqToTab(ofstream & outFile, FASTAinfo fastaseq){
	string line;
	line = fastaseq.seqHeader + "\t" + fastaseq.seq + "\t\n";
	outFile << line;
}


//scan through a fasta file or a tab delimited file and return the number of sequences and the lenght of the maximum sequence and length of the minimum esquence.
structvars::Fileinfo readfiles::readSeqLengths(const char* filename, string filetype,bool ignoreHeader){
	ifstream inFile;
	inFile.open(filename);
	structvars::Fileinfo seqdata;
	//structvars::FASTAinfo tempFastaSequence;
	structvars::FASTQinfo tempFastqSequence;
	int minLen = 1000, maxLen = 0,  totalLen = 0;
	int numSeqs = 0;
	string temp;

	if (filetype != "FASTA" && ignoreHeader)
		getline(inFile, temp);

	while (!inFile.eof()){
		
		if (filetype == "FASTA")
			tempFastqSequence = readIndividualSequence(inFile);
		else if(filetype=="FASTQ")
			tempFastqSequence = readIndividualFASTQSequence(inFile);
		else if (filetype == "TAB")
			tempFastqSequence = readIndividualTabSequence(inFile);
		else{
			printf("Incorrect filetype in readseqlengths function\n");
			exit(EXIT_FAILURE);
		}
		if (tempFastqSequence.seq != ""){//just make sure we actually read a proper line of sequence and sequence header
			numSeqs += 1;
			if (tempFastqSequence.seq.length() < minLen)
				minLen = tempFastqSequence.seq.length();
			if (tempFastqSequence.seq.length()>maxLen)
				maxLen = tempFastqSequence.seq.length();			
			AddSeqLenToVecDistribution(tempFastqSequence.seq.length(), seqdata.lengthDistribution);
			totalLen += tempFastqSequence.seq.length();
		}
		
	}
	printf("A total of %i sequences were found in the query file\n", numSeqs);
	inFile.close();
	seqdata.maxSeqLen = maxLen;
	seqdata.minSeqLen = minLen;
	seqdata.numSeqs = numSeqs;
	seqdata.avgSeqLen = (totalLen / numSeqs) + 1;
	sort(seqdata.lengthDistribution.begin(), seqdata.lengthDistribution.end(), [](const vector<int> & a, const vector<int> & b){ return (a[3] > b[3]); });

	return seqdata;
}

// trim whitespace from line
void readfiles::RemoveTrailingSpaces(string & s){
	s.erase(s.find_last_not_of(" \n\r\t") + 1);//remove from right side of string
	s.erase(0, s.find_first_not_of(" \n\r\t"));//remove from left side of stirng
}

//for this string:
//each variable is separated by commas:
//comma 1 = mindiag
//comma 2 = maxdiag
//comma 3 = fftMutations
//comma 4 = totalMutations
//comma 5 = numAlgn
//comma 6 = top_start
//comma 7 = bottom_start
//comma 8 = ins
//this variable defines WHERE there are insertion and deletions in the alignment.  Again this is another complex structure.  In this structure, each row is separated by ; and each column by %
//76%3;80%10 would be seen as:
//[76,3]
//[80,10]
//comma 9 = del
structvars::PairwiseClusterAlignment readfiles::ParsePairwiseAlignmentString(string pairwise_alignment_string){
	std::istringstream delim1, delim2, delim3; //delim1 will be for splitting string by commas, delim 2 will be for splitting string by ;, delim3 will be for splitting string by %
	structvars::PairwiseClusterAlignment tempvar;

	//first var
	delim1.clear();
	delim1.str(pairwise_alignment_string);
	string tempStrings, temprows, tempcols;
	vector<vector<int>> tempalignment;
	vector<int>temprowArray;
	int count_cols;
	int varsInArray[7];   //first seven elements will be easy to store in array
	int currentFieldCount = 0;
	for (int i = 0; i < 9; i++){
		getline(delim1, tempStrings, ','); //delimit string by commas
		if (i < 7){
			varsInArray[i] = std::stoi(tempStrings);
		}
		else{ //the last two variables are teh hardest to parse.  each semicolon is a row, and each % is a column
			delim2.clear();
			delim2.str(tempStrings);
			tempalignment.clear();
			while (getline(delim2, temprows, ';')){//delimit by semicolons (each of the insertion or deletion events in the alingment)
				delim3.clear();
				delim3.str(temprows);
				temprowArray.clear();
				count_cols = 0;
				while (getline(delim3, tempcols, '%')){
					temprowArray.push_back(std::stoi(tempcols));
					count_cols++;
				}
				temprowArray[0] -= 1; //matlab indexes start at 1 not 0, so adjut at this position
				tempalignment.push_back(temprowArray);
			}
			if (i == 7)
				tempvar.insertionEvents = tempalignment;
			else if (i == 8)
				tempvar.deletionEvents = tempalignment;
		}
	}
	//OK store the results into tempvar
	tempvar.minDiag = abs(varsInArray[0]);
	tempvar.maxDiag = abs(varsInArray[1]);
	tempvar.fftMutations = varsInArray[2];
	tempvar.totalMutations = varsInArray[3];
	//tempvar.numAlgn = varsInArray[4];
	//tempvar.top_start = varsInArray[5];
	//tempvar.bottom_start = varsInArray[6];
	return tempvar;
}

void readfiles::FindAbRegion(structvars::Abregion & ab_sub_region, string full_seq, int & startpos){
	std::size_t found;	
	if (ab_sub_region.sequence != ""){ //make sure there is even a sequence to be found
		found = full_seq.find(ab_sub_region.sequence, startpos);
		if (found != std::string::npos){			
			ab_sub_region.startpos = found;
			ab_sub_region.endpos = found + ab_sub_region.sequence.length() - 1;
			startpos = ab_sub_region.endpos + 1;
		}
		else{
			ab_sub_region.startpos = -1;
			ab_sub_region.endpos = -1;
			printf("There was a problem reading in the annotated subsequences");
			exit(EXIT_FAILURE);
		}
	}
	return;
}

