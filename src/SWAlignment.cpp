#include "SWAlignment.h"

using namespace std;
SWAlignment::SWAlignment(int max_seq_len, int max_germ_len){
	query = "";
	germline = "";
	maxSeqInsSize = max_seq_len + max_germ_len;// 5 * swParams.maxGap;
	maxSeqDelSize = max_seq_len + max_germ_len;// 5 * swParams.maxGap;
	maxScore = 0;
	gapPenalties = 0;
	seqLen = max_seq_len;
	germLen = max_germ_len;
	seqDel = new int*[maxSeqDelSize];
	seqIns = new int*[maxSeqInsSize];
	for (int k = 0; k < maxSeqDelSize; k++){
		seqDel[k] = new int[2];
	}
	for (int k = 0; k < maxSeqInsSize; k++){
		seqIns[k] = new int[2];
	}

	alignmentPath = (int ***)malloc(3 * sizeof(int**)); //make alignment path for moving 1) DIAG, 2) RIGHT, AND 3) DOWN
	scoreMatrix = (double ***)malloc(3 * sizeof(double**));

	for (int j = 0; j < 3; j++){
		scoreMatrix[j] = (double**)malloc((seqLen + 2)* sizeof(double*));
		alignmentPath[j] = (int**)malloc((seqLen + 2)*sizeof(int*));
		for (int i = 0; i < seqLen + 2; i++){
			scoreMatrix[j][i] = (double*)calloc((germLen + 1), sizeof(double));
			alignmentPath[j][i] = (int*)calloc((germLen + 1), sizeof(int));
		}
	}

	tracebackData = (double**)calloc((germLen + 1), sizeof(double*));
	for (int j = 0; j < germLen + 1; j++){
		tracebackData[j] = (double*)calloc(7, sizeof(double));
	}

	int maxPossibleLen = max(seqLen, germLen) + 1;
	gaplessScores = new double[maxPossibleLen];

	final_posIJ[0] = 0;
	final_posIJ[1] = 0;
	alignedGerm.resize(germLen + seqLen);
	alignedQuery.resize(germLen + seqLen);
	gaplessGermline.resize(germLen + seqLen);
	gaplessQuery.resize(germLen + seqLen);
}

SWAlignment::SWAlignment(int max_seq_len, int max_germ_len, structvars::SWAlignSettings settings){
	swParams = settings;
	query = "";
	germline = "";

	maxSeqInsSize = 5 * swParams.maxGap;
	maxSeqDelSize = 5 * swParams.maxGap;
	maxScore = 0;
	gapPenalties = 0;
	seqLen = max_seq_len;
	germLen = max_germ_len;
	seqDel = new int*[maxSeqDelSize];
	seqIns = new int*[maxSeqInsSize];
	for (int k = 0; k < maxSeqDelSize; k++){
		seqDel[k] = new int[2];
	}
	for (int k = 0; k < maxSeqInsSize; k++){
		seqIns[k] = new int[2];
	}

	alignmentPath = (int ***)malloc(3 * sizeof(int**)); //make alignment path for moving 1) DIAG, 2) RIGHT, AND 3) DOWN
	scoreMatrix = (double ***)malloc(3 * sizeof(double**));

	for (int j = 0; j < 3; j++){
		scoreMatrix[j] = (double**)malloc((seqLen + 2)* sizeof(double*));
		alignmentPath[j] = (int**)malloc((seqLen + 2)*sizeof(int*));
		for (int i = 0; i < seqLen + 2; i++){
			scoreMatrix[j][i] = (double*)calloc((germLen + 1), sizeof(double));
			alignmentPath[j][i] = (int*)calloc((germLen + 1), sizeof(int));
		}

	}

	tracebackData = (double**)calloc((germLen + 1), sizeof(double*));
	for (int j = 0; j < germLen + 1; j++){
		tracebackData[j] = (double*)calloc(7, sizeof(double));
	}

	int maxPossibleLen = max(seqLen, germLen) + 1;
	gaplessScores = new double[maxPossibleLen];

	final_posIJ[0] = 0;
	final_posIJ[1] = 0;
	alignedGerm.resize(germLen + seqLen);
	alignedQuery.resize(germLen + seqLen);
	gaplessGermline.resize(germLen + seqLen);
	gaplessQuery.resize(germLen + seqLen);
}

SWAlignment::~SWAlignment(){
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < seqLen + 2; j++){
			free(scoreMatrix[i][j]);
			free(alignmentPath[i][j]);
		}
		free(scoreMatrix[i]);
		free(alignmentPath[i]);
	}
	free(scoreMatrix);
	free(alignmentPath);

	for (int i = 0; i < germLen + 1; i++){
		free(tracebackData[i]);
	}
	free(tracebackData);

	for (int i = 0; i < maxSeqDelSize; i++)
		delete[] seqDel[i];
	for (int i = 0; i < maxSeqInsSize; i++)
		delete[] seqIns[i];

	delete[] seqDel;
	delete[] seqIns;
	delete[] gaplessScores;
}


//This function will perform a gapless alignment (not allow any other paths other than the path defined by the diagonal in the variable info
//info[0] = > seqStart
//info[1] = > germStart, 
//info[2] = > alignment length
void SWAlignment::GaplessAlignDiag(const std::string & sequence, const std::string & germline, int *info){
	//int lenCompare = (info[0] + info[2] - 1 > sequence.length() - 1 || info[1] + info[2] - 1 > germline.length() - 1) ? std::min(sequence.length() - info[0], germline.length() - info[1]) : info[2];// min(sequence.length() - info[0] , germline.length() - (info[1] + info[2]) - 1);	
	int lenCompare = std::min(sequence.length() - info[0], germline.length() - info[1]);// : info[2];// min(sequence.length() - info[0] , germline.length() - (info[1] + info[2]) - 1);
	lenCompare = std::min(info[2], lenCompare);

	int score = 0;
	maxScore = 0;
	gaplessScores[0] = 0;

	int end = lenCompare, start = 0;

	for (int i = 0; i < lenCompare; i++){
		score = (sequence[i + info[0]] == NOT_NT || germline[i + info[1]] == NOT_NT) ? 0 : ((sequence[i + info[0]] == germline[i + info[1]]) ? swParams.matchScore : -1 * abs(swParams.mismatchScore));
		score = score + gaplessScores[i];
		gaplessScores[i + 1] = score <= 0 ? 0 : score;

		if (score > maxScore){
			maxScore = score;
			end = i + 1;
		}
	}
	start = end;

	while (gaplessScores[start] > 0)
			start--;
	start += 1;

	info[2] = end - start;
	info[0] += start;
	info[1] += start;

	seqEnd = info[0] + info[2] - 1;
	germEnd = info[1] + info[2] - 1;
	seqStart = info[0];
	germStart = info[1];
}


//THIS FUNCTION WILL ALIGN A QUERY SEQUENCE TO THE GERMLINE SEQUENCE AS DEFINED BY THE GERMLINE COORINATES.
//So the germline sequence position will be the x-axis (columns) and query sequence will be the y-axis (rows)
//coordinates = > [0] => germ start, [1] => germ End, [2] => diagonal, [3] = > min diag, [4] = > max diag
//scores => [0] = match score, [1] = mismatch penalty, [2] = gap open penalty, [3] = gap extend penalty
//query nucleotide sequence
//germline nucleotide sequence
//the diagonal is defined as sequence position - germline position
//"mindiagonal" is defined with respect to sequence and "max diagonal" is defined with respect to sequence
void SWAlignment::SWAlignDiag(const string & seq, const string & germSeq, int *coordinates) //, double posIJ[])
{
	int minDiag = coordinates[3], maxDiag = coordinates[4], bestDiag = coordinates[2], germStart = coordinates[0], germEnd = coordinates[1];
	query = seq;
	germline = germSeq;
	double match = swParams.matchScore, mismatch = abs(swParams.mismatchScore), open_gap = swParams.swGapOpen, extend_gap = swParams.swExtendGap;
	double neg_inf = std::numeric_limits<double>::lowest();

	int g = coordinates[6];// germSeq.length();
	int s = seq.length();
	int seqStart = coordinates[5];// 0;

	//germStart = std::max(germStart - maxDiag - 1, 0) + 1; //maxdiag again is with respect to sequence not germline, so moving horinzontally left can only go as far as max diag
	//germEnd = std::min(germEnd + minDiag + 1, g - 1) + 1; //see above
	germStart = germStart + 1;
	germEnd = germEnd + 1;

	int startingDiag = 0;
	double currentScore = 0.0, current_best_score = 0.0, match_score;

	int currentInd;

	bool is_equal, is_nt;

	double each_dir_score[] = { 0, 0, 0, 0 }; //each index corresponds to movements in alignment path (SEE ENUM DECLARED ABOVE to know which movements correspond to which direction)

	double *best_score, max_alignment_score = 0;
	int best_dir;

	int starting_query_index, min_query_index, max_query_index;
	double penalty_horizontal, penalty_diagonal, penalty_vertical;
	int seqGap, germGap;

	//initialize the matrix values upstream of the start position to 0
	min_query_index = std::min(s-1,std::max(seqStart, (((germStart - 1) + bestDiag) - minDiag))) + 1;
	max_query_index = std::max(seqStart,std::min(s - 1, (((germStart - 1) + bestDiag) + maxDiag))) + 1;

	scoreMatrix[DIAG][min_query_index - 1][germStart - 1] = 0;
	scoreMatrix[RIGHT][min_query_index - 1][germStart - 1] = 0;
	scoreMatrix[DOWN][min_query_index - 1][germStart - 1] = 0;

	alignmentPath[DIAG][min_query_index - 1][germStart - 1] = START;
	alignmentPath[RIGHT][min_query_index - 1][germStart - 1] = START;
	alignmentPath[DOWN][min_query_index - 1][germStart - 1] = START;

	scoreMatrix[DIAG][max_query_index + 1][germStart - 1] = 0;
	scoreMatrix[RIGHT][max_query_index + 1][germStart - 1] = neg_inf; //make this event/occurence/movement impossible
	alignmentPath[DIAG][max_query_index + 1][germStart - 1] = START;
	alignmentPath[RIGHT][max_query_index + 1][germStart - 1] = START;
	alignmentPath[DOWN][max_query_index + 1][germStart - 1] = START;

	for (int i = min_query_index; i <= max_query_index; i++){
		scoreMatrix[DIAG][i][germStart - 1] = 0;
		scoreMatrix[RIGHT][i][germStart - 1] = neg_inf; //make this event/occurence/movement impossible

		alignmentPath[DIAG][i][germStart - 1] = START;
		alignmentPath[RIGHT][i][germStart - 1] = START;
		alignmentPath[DOWN][i][germStart - 1] = START;
	}

	final_posIJ[0] = 0;
	final_posIJ[1] = 0;

	/*INITALIZING DONE*/

	for (int j = germStart; j <= germEnd; j++){
		starting_query_index = j + bestDiag; //sequence = germline+diagonal

		//initialize the matrix values upstream of the start position to 0
		min_query_index = std::min(s-1,std::max(seqStart, (starting_query_index - 1 - minDiag))) + 1;
		max_query_index = std::max(seqStart,std::min(s - 1, (starting_query_index - 1 + maxDiag))) + 1;

		if (min_query_index <= max_query_index){
			scoreMatrix[DIAG][min_query_index - 1][j] = 0;
			scoreMatrix[DOWN][min_query_index - 1][j] = neg_inf;

			alignmentPath[DIAG][min_query_index - 1][j] = START;
			alignmentPath[DOWN][min_query_index - 1][j] = START;
			alignmentPath[RIGHT][min_query_index - 1][j] = START;

			scoreMatrix[DIAG][max_query_index + 1][j] = 0;
			scoreMatrix[RIGHT][max_query_index + 1][j] = neg_inf;

			alignmentPath[DIAG][max_query_index + 1][j] = START; //pad row right below last possible diagonal with 0 and 3 respectively
			alignmentPath[DOWN][max_query_index + 1][j] = START; //pad row right below last possible diagonal with 0 and 3 respectively
			alignmentPath[RIGHT][max_query_index + 1][j] = START; //pad row right below last possible diagonal with 0 and 3 respectively
		}
		for (int i = min_query_index; i <= max_query_index; i++)
		{
			//the following values correspond to a specific movement: {0,1,2} corresponds to {moving right, moving down, and moving diagonally} along matrix

			//key thing to note: first row and column of score matrix are "0" or "NULL", so this meanst that row[1] of score matrix will correspond to char[0] of query and column [1] of score matrix will correspond to char[0] of germline

			each_dir_score[START] = 0.0;
			current_best_score = 0.0;
			best_dir = 0;
			is_nt = (seq[i - 1] != NOT_NT && germSeq[j - 1] != NOT_NT); //checks whether both nucleotides are .. nucleotides (not not unkonw characters i.e. N or X

			match_score = seq[i - 1] == germSeq[j - 1] ? is_nt*match : -is_nt*mismatch;

			penalty_diagonal = scoreMatrix[DIAG][i][j - 1] - open_gap;
			penalty_horizontal = scoreMatrix[RIGHT][i][j - 1] - extend_gap;
			if (penalty_diagonal > penalty_horizontal){
				scoreMatrix[RIGHT][i][j] = penalty_diagonal;
				alignmentPath[RIGHT][i][j] = DIAG;
			}
			else{
				scoreMatrix[RIGHT][i][j] = penalty_horizontal;
				alignmentPath[RIGHT][i][j] = RIGHT;
			}

			penalty_diagonal = scoreMatrix[DIAG][i - 1][j] - open_gap;
			penalty_vertical = scoreMatrix[DOWN][i - 1][j] - extend_gap;

			if (penalty_diagonal > penalty_vertical){
				scoreMatrix[DOWN][i][j] = penalty_diagonal;
				alignmentPath[DOWN][i][j] = DIAG;
			}
			else{
				scoreMatrix[DOWN][i][j] = penalty_vertical;
				alignmentPath[DOWN][i][j] = DOWN;
			}

			best_dir = scoreMatrix[RIGHT][i - 1][j - 1] > scoreMatrix[DIAG][i - 1][j - 1] ? RIGHT : DIAG;

			best_dir = scoreMatrix[DOWN][i - 1][j - 1] > scoreMatrix[best_dir][i - 1][j - 1] ? DOWN : best_dir;

			current_best_score = scoreMatrix[best_dir][i - 1][j - 1] + match_score;
			if (0 >= current_best_score){
				best_dir = START;
				current_best_score = 0;
			}

			scoreMatrix[DIAG][i][j] = current_best_score;//*best_score;
			alignmentPath[DIAG][i][j] = best_dir;

			if (current_best_score > max_alignment_score){// *best_score //always expect the final best score to end on the diagonal (no gaps at teh end that is)
				final_posIJ[0] = i;
				final_posIJ[1] = j;
				max_alignment_score = current_best_score;// *best_score;
			}
		}
	}
	maxScore = max_alignment_score;
	return;
}


//Function for overlap alignment
//THIS FUNCTION WILL ALIGN A QUERY SEQUENCE TO THE GERMLINE SEQUENCE AS DEFINED BY THE GERMLINE COORINATES.
//So the germline sequence position will be the x-axis (columns) and query sequence will be the y-axis (rows)
//coordinates = > [0] => germ start, [1] => germ End, [2] => diagonal, [3] = > min diag, [4] = > max diag, [5]=> queryStart, [6] = > germ end
//scores => [0] = match score, [1] = mismatch penalty, [2] = gap open penalty, [3] = gap extend penalty
//query nucleotide sequence
//germline nucleotide sequence
//the diagonal is defined as sequence position - germline position
//"mindiagonal" is defined with respect to sequence and "max diagonal" is defined with respect to sequence
void SWAlignment::OverlapAlignDiag(const string & seq, const string & germSeq, int *coordinates) //, double posIJ[])
{
	int minDiag = coordinates[3], maxDiag = coordinates[4], bestDiag = coordinates[2], germStart = coordinates[0], germEnd = coordinates[1];
	query = seq;
	germline = germSeq;
	double match = swParams.matchScore, mismatch = abs(swParams.mismatchScore), open_gap = swParams.swGapOpen, extend_gap = swParams.swExtendGap;
	double neg_inf = std::numeric_limits<double>::lowest();
	int g = coordinates[6];// germSeq.length();
	int s = seq.length();
	int seqStart = coordinates[5];// 0;

	
	//germStart = std::max(germStart - maxDiag - 1, 0) + 1; //maxdiag again is with respect to sequence not germline, so moving horinzontally left can only go as far as max diag 
	germStart = germStart + 1; //line above commented out....
	germEnd = germEnd + 1;// std::min(germEnd + minDiag + 1, g - 1) + 1; //see above

	int startingDiag = 0;
	double currentScore = 0.0, current_best_score = 0.0, match_score;

	int currentInd;

	bool is_equal, is_nt;

	double each_dir_score[] = { 0, 0, 0, 0 }; //each index corresponds to movements in alignment path (SEE ENUM DECLARED ABOVE to know which movements correspond to which direction)

	double *best_score, max_alignment_score = neg_inf;
	int best_dir;

	int starting_query_index, min_query_index, max_query_index;
	double penalty_horizontal, penalty_diagonal, penalty_vertical;
	int seqGap, germGap;

	//initialize the matrix values upstream of the start position to 0
	min_query_index = std::min(s-1,std::max(seqStart, (((germStart - 1) + bestDiag) - minDiag)))+1 ;
	//min_query_index = std::min(min_query_index, s - 1);

	max_query_index = std::max(seqStart,std::min(s - 1, (((germStart - 1) + bestDiag) + maxDiag)))+1;
	//max_query_index = std::max(seqStart, max_query_index);

	//min_query_index += 1;
	//max_query_index += 1;

	for (int i = min_query_index - 1; i <= max_query_index + 1; i++){
		scoreMatrix[DIAG][i][germStart - 1] = 0;
		scoreMatrix[RIGHT][i][germStart - 1] = neg_inf; //make this event/occurence/movement impossible
		scoreMatrix[DOWN][i][germStart - 1] = neg_inf; //make this event/occurence/movement impossible

		alignmentPath[DIAG][i][germStart - 1] = START;
		alignmentPath[RIGHT][i][germStart - 1] = START;
		alignmentPath[DOWN][i][germStart - 1] = START;
	}

	final_posIJ[0] = 0;
	final_posIJ[1] = 0;

	/*INITALIZING DONE*/

	for (int j = germStart; j <= germEnd; j++){
		starting_query_index = j + bestDiag; //sequence = germline+diagonal

		//initialize the matrix values upstream of the start position to 0
		min_query_index = std::min(s-1,std::max(seqStart, (starting_query_index - 1 - minDiag))) + 1;
		max_query_index = std::max(seqStart,std::min(s - 1, (starting_query_index - 1 + maxDiag))) + 1;

		if (min_query_index <= max_query_index){
			scoreMatrix[DIAG][min_query_index - 1][j] = 0;
			scoreMatrix[DOWN][min_query_index - 1][j] = neg_inf;//std::numeric_limits<double>::lowest();

			alignmentPath[DIAG][min_query_index - 1][j] = START;
			alignmentPath[DOWN][min_query_index - 1][j] = START;
			alignmentPath[RIGHT][min_query_index - 1][j] = START;

			scoreMatrix[DIAG][max_query_index + 1][j] = 0;
			scoreMatrix[RIGHT][max_query_index + 1][j] = neg_inf;// std::numeric_limits<double>::lowest();

			alignmentPath[DIAG][max_query_index + 1][j] = START; //pad row right below last possible diagonal with 0 and 3 respectively
			alignmentPath[DOWN][max_query_index + 1][j] = START; //pad row right below last possible diagonal with 0 and 3 respectively
			alignmentPath[RIGHT][max_query_index + 1][j] = START; //pad row right below last possible diagonal with 0 and 3 respectively
		}
		for (int i = min_query_index; i <= max_query_index; i++)
		{
			//the following values correspond to a specific movement: {0,1,2} corresponds to {moving right, moving down, and moving diagonally} along matrix

			//key thing to note: first row and column of score matrix are "0" or "NULL", so this meanst that row[1] of score matrix will correspond to char[0] of query and column [1] of score matrix will correspond to char[0] of germline

			each_dir_score[START] = 0.0;
			current_best_score = 0.0;
			best_dir = 0;
			is_nt = (seq[i - 1] != NOT_NT && germSeq[j - 1] != NOT_NT); //checks whether both nucleotides are .. nucleotides (not not unkonw characters i.e. N or X

			match_score = seq[i - 1] == germSeq[j - 1] ? is_nt*match : -is_nt*mismatch;

			penalty_diagonal = scoreMatrix[DIAG][i][j - 1] - open_gap;
			penalty_horizontal = scoreMatrix[RIGHT][i][j - 1] - extend_gap;
			if (penalty_diagonal > penalty_horizontal){
				scoreMatrix[RIGHT][i][j] = penalty_diagonal;
				alignmentPath[RIGHT][i][j] = DIAG;
			}
			else{
				scoreMatrix[RIGHT][i][j] = penalty_horizontal;
				alignmentPath[RIGHT][i][j] = RIGHT;
			}

			penalty_diagonal = scoreMatrix[DIAG][i - 1][j] - open_gap;
			penalty_vertical = scoreMatrix[DOWN][i - 1][j] - extend_gap;

			if (penalty_diagonal > penalty_vertical){
				scoreMatrix[DOWN][i][j] = penalty_diagonal;
				alignmentPath[DOWN][i][j] = DIAG;
			}
			else{
				scoreMatrix[DOWN][i][j] = penalty_vertical;
				alignmentPath[DOWN][i][j] = DOWN;
			}

			best_dir = scoreMatrix[RIGHT][i - 1][j - 1] > scoreMatrix[DIAG][i - 1][j - 1] ? RIGHT : DIAG;

			best_dir = scoreMatrix[DOWN][i - 1][j - 1] > scoreMatrix[best_dir][i - 1][j - 1] ? DOWN : best_dir;

			current_best_score = scoreMatrix[best_dir][i - 1][j - 1] + match_score;
			//if (0 >= current_best_score){
			//	best_dir = START;
			//	current_best_score = 0;
			//}

			scoreMatrix[DIAG][i][j] = current_best_score;//*best_score;
			alignmentPath[DIAG][i][j] = best_dir;

			if ((j == germEnd || i == s) && (current_best_score > max_alignment_score)){//using overlap alignment rules, so find best alignment at the end/edges of the alignment matrix
				final_posIJ[0] = i;
				final_posIJ[1] = j;
				max_alignment_score = current_best_score;// *best_score;
			}
		}
	}
	maxScore = max_alignment_score;
	//Traceback(alignmentPath);
	return;
}



//performs an entire alignment using both sequences (does not make any assumptions about min or max diags)
//this will be a slower function as opposed to functions such as OverlapAlignDiag where the variables and sizes are predefined
void SWAlignment::OverlapAlignComplete(const std::string & seq1, const std::string & seqGerm){
	int numXCoord = seq1.length();
	int numYCoord = seqGerm.length();
	int ***tempAlignmentPath;
	double ***tempScoreMatrix;
	double neg_inf = std::numeric_limits<double>::lowest();
	tempAlignmentPath = (int ***)malloc(3 * sizeof(int**)); //make alignment path for moving 1) DIAG, 2) RIGHT, AND 3) DOWN
	tempScoreMatrix = (double ***)malloc(3 * sizeof(double**));
	query = seq1;
	germline = seqGerm;
	germLen = numYCoord;
	seqLen = numXCoord;

	for (int j = 0; j < 3; j++){
		tempScoreMatrix[j] = (double**)malloc((numXCoord + 1)* sizeof(double*));
		tempAlignmentPath[j] = (int**)malloc((numXCoord + 1)*sizeof(int*));
		for (int i = 0; i < numXCoord + 1; i++){
			tempScoreMatrix[j][i] = (double*)malloc((numYCoord + 1)*sizeof(double));
			tempAlignmentPath[j][i] = (int*)malloc((numYCoord + 1)*sizeof(int));
		}
	}
	double each_dir_score[] = { 0, 0, 0, 0 }; //each index corresponds to movements in alignment path (SEE ENUM DECLARED ABOVE to know which movements correspond to which direction)
	double current_best_score = 0.0, currentScore = 0.0, match_score;

	int startingDiag = 0, currentInd, best_dir;

	bool is_equal, is_nt;

	double *best_score, max_alignment_score = neg_inf;

	double penalty_horizontal, penalty_diagonal, penalty_vertical;
	int seqGap, germGap;
	double match = swParams.matchScore, mismatch = abs(swParams.mismatchScore), open_gap = swParams.swGapOpen, extend_gap = swParams.swExtendGap;

	//intialize the first row and column: 
	for (int i = 0; i < numXCoord + 1; i++){
		tempScoreMatrix[DIAG][i][0] = 0;
		tempScoreMatrix[RIGHT][i][0] = neg_inf;
		tempScoreMatrix[DOWN][i][0] = neg_inf;
		tempAlignmentPath[DIAG][i][0] = START;
		tempAlignmentPath[RIGHT][i][0] = START;
		tempAlignmentPath[DOWN][i][0] = START;
	}

	for (int j = 0; j < numYCoord + 1; j++){
		tempScoreMatrix[DIAG][0][j] = 0;
		tempScoreMatrix[RIGHT][0][j] = neg_inf;
		tempScoreMatrix[DOWN][0][j] = neg_inf;
		tempAlignmentPath[DIAG][0][j] = START;
		tempAlignmentPath[RIGHT][0][j] = START;
		tempAlignmentPath[DOWN][0][j] = START;
	}
	/*initializing done*/

	final_posIJ[0] = 0;
	final_posIJ[1] = 0;


	for (int i = 1; i < numXCoord + 1; i++){
		for (int j = 1; j < numYCoord + 1; j++){
			each_dir_score[START] = 0.0;
			current_best_score = 0.0;
			best_dir = 0;
			is_nt = (seq1[i - 1] != NOT_NT && seqGerm[j - 1] != NOT_NT); //checks whether both nucleotides are .. nucleotides (not not unkonw characters i.e. N or X

			match_score = seq1[i - 1] == seqGerm[j - 1] ? is_nt*match : -is_nt*mismatch;

			penalty_diagonal = tempScoreMatrix[DIAG][i][j - 1] - open_gap;
			penalty_horizontal = tempScoreMatrix[RIGHT][i][j - 1] - extend_gap;
			if (penalty_diagonal > penalty_horizontal){
				tempScoreMatrix[RIGHT][i][j] = penalty_diagonal;
				tempAlignmentPath[RIGHT][i][j] = DIAG;
			}
			else{
				tempScoreMatrix[RIGHT][i][j] = penalty_horizontal;
				tempAlignmentPath[RIGHT][i][j] = RIGHT;
			}

			penalty_diagonal = tempScoreMatrix[DIAG][i - 1][j] - open_gap;
			penalty_vertical = tempScoreMatrix[DOWN][i - 1][j] - extend_gap;

			if (penalty_diagonal > penalty_vertical){
				tempScoreMatrix[DOWN][i][j] = penalty_diagonal;
				tempAlignmentPath[DOWN][i][j] = DIAG;
			}
			else{
				tempScoreMatrix[DOWN][i][j] = penalty_vertical;
				tempAlignmentPath[DOWN][i][j] = DOWN;
			}

			best_dir = tempScoreMatrix[RIGHT][i - 1][j - 1] > tempScoreMatrix[DIAG][i - 1][j - 1] ? RIGHT : DIAG;

			best_dir = tempScoreMatrix[DOWN][i - 1][j - 1] > tempScoreMatrix[best_dir][i - 1][j - 1] ? DOWN : best_dir;

			current_best_score = tempScoreMatrix[best_dir][i - 1][j - 1] + match_score;
			//if (0 >= current_best_score){
			//	best_dir = START;
			//	current_best_score = 0;
			//}

			tempScoreMatrix[DIAG][i][j] = current_best_score;//*best_score;
			tempAlignmentPath[DIAG][i][j] = best_dir;

			if ((j == numYCoord || i == numXCoord) && (current_best_score > max_alignment_score)){//using overlap alignment rules, so find best alignment at the end/edges of the alignment matrix
				final_posIJ[0] = i;
				final_posIJ[1] = j;
				max_alignment_score = current_best_score;// *best_score;
			}
		}
	}

	maxScore = max_alignment_score;

	Traceback(tempAlignmentPath);

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < numXCoord + 1; j++){
			free(tempScoreMatrix[i][j]);
			free(tempAlignmentPath[i][j]);
		}
		free(tempScoreMatrix[i]);
		free(tempAlignmentPath[i]);
	}
	free(tempScoreMatrix);
	free(tempAlignmentPath);



	return;
}


void SWAlignment::Traceback(int ***extAlgnPath){
	string endSeq, begSeq;
	alignedGerm = "";
	alignedQuery = "";
	gaplessGermline = "";
	gaplessQuery = "";
	gapPenalties = 0;
	bool delOpen = false, insOpen = false;

	int k = germLen + seqLen - 1;
	int s = k, g = k;

	alignedGerm.resize(k + 1);
	alignedQuery.resize(k + 1);
	gaplessGermline.resize(k + 1);
	gaplessQuery.resize(k + 1);

	int seq_pos = final_posIJ[0];
	int germ_pos = final_posIJ[1];
	int next_matrix = DIAG, current_matrix, prev_matrix;//always start at the diag position
	seqEnd = seq_pos - 1;
	germEnd = germ_pos - 1;
	numSeqIns[0] = 0;
	numSeqDel[0] = 0;
	numSeqDel[1] = 0;
	numSeqIns[1] = 0;
	int **tempSeqIns, **tempSeqDel;

	endSeq = query.substr(seq_pos, query.length() - seq_pos);
	bool beg_of_seq = false;

	do{
		current_matrix = next_matrix;

		next_matrix = extAlgnPath[current_matrix][seq_pos][germ_pos];

		if (next_matrix == START)
			break;
		switch (current_matrix){
		case DIAG:
			alignedGerm[k] = germline[germ_pos - 1];
			alignedQuery[k] = query[seq_pos - 1];
			gaplessGermline[g] = germline[germ_pos - 1];
			gaplessQuery[s] = query[seq_pos - 1];
			s--;
			g--;
			seq_pos--;
			germ_pos--;
			break;

		case RIGHT:
			alignedGerm[k] = germline[germ_pos - 1];
			alignedQuery[k] = '-';
			gaplessQuery[s] = NOT_NT; //add an X in the seuence, but do not add any letter ot aligned germline

			numSeqIns[0] += 1;
			if (prev_matrix == RIGHT){
				gapPenalties += swParams.swExtendGap;
				seqIns[numSeqIns[1] - 1][1] += 1; //add 1 to the "seq insertions" list
			}
			else{
				gapPenalties += swParams.swGapOpen;
				seqIns[numSeqIns[1]][1] = 1; //add 1 to the "seq insertions" list
				numSeqIns[1]++;
			}
			seqIns[numSeqIns[1] - 1][0] = seq_pos - 1;
			germ_pos--;
			s--;
			break;

		case DOWN:
			alignedGerm[k] = '-';
			gaplessGermline[g] = NOT_NT;
			alignedQuery[k] = query[seq_pos - 1];

			numSeqDel[0] += 1;
			if (prev_matrix == DOWN){
				gapPenalties += swParams.swExtendGap;
				seqDel[numSeqDel[1] - 1][1] += 1;
			}
			else{
				gapPenalties += swParams.swGapOpen;
				seqDel[numSeqDel[1]][1] = 1;
				numSeqDel[1]++;
			}
			seqDel[numSeqDel[1] - 1][0] = seq_pos - 1;
			seq_pos--;
			g--;
			break;
		}
		prev_matrix = current_matrix;
		k--;
	} while (next_matrix != START);
	alignedGerm.erase(0, k + 1);
	alignedQuery.erase(0, k + 1);

	seqStart = seq_pos;
	germStart = germ_pos;

	if (alignedQuery[0] == '-'){
		numSeqIns[0] -= seqIns[numSeqIns[1] - 1][1];
		double scoreAdjust = swParams.swExtendGap*(seqIns[numSeqIns[1] - 1][1] - 1) + swParams.swGapOpen;
		gapPenalties -= scoreAdjust;
		maxScore += scoreAdjust;
		germStart += seqIns[numSeqIns[1] - 1][1];
		s += seqIns[numSeqIns[1] - 1][1];
		alignedQuery.erase(0, seqIns[numSeqIns[1] - 1][1]);
		alignedGerm.erase(0, seqIns[numSeqIns[1] - 1][1]);
		seqIns[numSeqIns[1] - 1][1] = 0;
		numSeqIns[1] -= 1;
	}
	if (alignedGerm[0] == '-'){
		numSeqDel[0] -= seqDel[numSeqDel[1] - 1][1];
		double scoreAdjust = swParams.swExtendGap*(seqDel[numSeqDel[1] - 1][1] - 1) + swParams.swGapOpen;
		gapPenalties -= scoreAdjust;
		maxScore += scoreAdjust;
		seqStart += seqDel[numSeqDel[1] - 1][1];
		g += seqDel[numSeqDel[1] - 1][1];
		alignedQuery.erase(0, seqDel[numSeqDel[1] - 1][1]);
		alignedGerm.erase(0, seqDel[numSeqDel[1] - 1][1]);

		seqDel[numSeqDel[1] - 1][1] = 0;
		numSeqDel[1] -= 1;
	}

	begSeq = query.substr(0, seqStart);

	gaplessGermline.erase(0, g + 1);//this could be wrong
	gaplessQuery.erase(0, s + 1);
	gaplessQuery = begSeq + gaplessQuery + endSeq;
}


void SWAlignment::Traceback(){
	string endSeq, begSeq;
	alignedGerm = "";
	alignedQuery = "";
	gaplessGermline = "";
	gaplessQuery = "";
	gapPenalties = 0;
	bool delOpen = false, insOpen = false;

	int k = germLen + seqLen - 1;
	int s = k, g = k;

	alignedGerm.resize(k + 1);
	alignedQuery.resize(k + 1);
	gaplessGermline.resize(k + 1);
	gaplessQuery.resize(k + 1);

	int seq_pos = final_posIJ[0];
	int germ_pos = final_posIJ[1];
	int next_matrix = DIAG, current_matrix, prev_matrix;//always start at the diag position
	seqEnd = seq_pos - 1;
	germEnd = germ_pos - 1;
	numSeqIns[0] = 0;
	numSeqDel[0] = 0;
	numSeqDel[1] = 0;
	numSeqIns[1] = 0;

	endSeq = query.substr(seq_pos, query.length() - seq_pos);
	bool beg_of_seq = false;

	do{
		current_matrix = next_matrix;

		next_matrix = alignmentPath[current_matrix][seq_pos][germ_pos];

		if (next_matrix == START)
			break;

		switch (current_matrix){
		case DIAG:
			alignedGerm[k] = germline[germ_pos - 1];
			alignedQuery[k] = query[seq_pos - 1];
			gaplessGermline[g] = germline[germ_pos - 1];
			gaplessQuery[s] = query[seq_pos - 1];
			s--;
			g--;
			seq_pos--;
			germ_pos--;
			break;

		case RIGHT:
			alignedGerm[k] = germline[germ_pos - 1];
			alignedQuery[k] = '-';
			gaplessQuery[s] = NOT_NT; //add an X in the seuence, but do not add any letter ot aligned germline

			numSeqIns[0] += 1;
			if (prev_matrix == RIGHT){
				gapPenalties += swParams.swExtendGap;
				seqIns[numSeqIns[1] - 1][1] += 1; //add 1 to the "seq insertions" list
			}
			else{
				gapPenalties += swParams.swGapOpen;
				seqIns[numSeqIns[1]][1] = 1; //add 1 to the "seq insertions" list
				numSeqIns[1]++;
			}
			seqIns[numSeqIns[1] - 1][0] = seq_pos - 1;
			germ_pos--;
			s--;
			break;

		case DOWN:
			alignedGerm[k] = '-';
			gaplessGermline[g] = NOT_NT;
			alignedQuery[k] = query[seq_pos - 1];

			numSeqDel[0] += 1;
			if (prev_matrix == DOWN){
				gapPenalties += swParams.swExtendGap;
				seqDel[numSeqDel[1] - 1][1] += 1;
			}
			else{
				gapPenalties += swParams.swGapOpen;
				seqDel[numSeqDel[1]][1] = 1;
				numSeqDel[1]++;
			}
			seqDel[numSeqDel[1] - 1][0] = seq_pos - 1;
			seq_pos--;
			g--;
			break;
		}
		prev_matrix = current_matrix;
		k--;
	} while (next_matrix != START);
	alignedGerm.erase(0, k + 1);
	alignedQuery.erase(0, k + 1);

	seqStart = seq_pos;
	germStart = germ_pos;

	if (alignedQuery[0] == '-'){
		numSeqIns[0] -= seqIns[numSeqIns[1] - 1][1];
		double scoreAdjust = swParams.swExtendGap*(seqIns[numSeqIns[1] - 1][1] - 1) + swParams.swGapOpen;
		gapPenalties -= scoreAdjust;
		maxScore += scoreAdjust;
		germStart += seqIns[numSeqIns[1] - 1][1];
		s += seqIns[numSeqIns[1] - 1][1];
		alignedQuery.erase(0, seqIns[numSeqIns[1] - 1][1]);
		alignedGerm.erase(0, seqIns[numSeqIns[1] - 1][1]);
		seqIns[numSeqIns[1] - 1][1] = 0;
		numSeqIns[1] -= 1;
	}
	if (alignedGerm[0] == '-'){
		numSeqDel[0] -= seqDel[numSeqDel[1] - 1][1];
		double scoreAdjust = swParams.swExtendGap*(seqDel[numSeqDel[1] - 1][1] - 1) + swParams.swGapOpen;
		gapPenalties -= scoreAdjust;
		maxScore += scoreAdjust;
		seqStart += seqDel[numSeqDel[1] - 1][1];
		g += seqDel[numSeqDel[1] - 1][1];
		alignedQuery.erase(0, seqDel[numSeqDel[1] - 1][1]);
		alignedGerm.erase(0, seqDel[numSeqDel[1] - 1][1]);

		seqDel[numSeqDel[1] - 1][1] = 0;
		numSeqDel[1] -= 1;
	}

	begSeq = query.substr(0, seqStart);

	gaplessGermline.erase(0, g + 1);
	gaplessQuery.erase(0, s + 1);
	gaplessQuery = begSeq + gaplessQuery + endSeq;
}


void SWAlignment::TracebackWithAnnotation(){
	string endSeq, begSeq;
	alignedGerm = "";
	alignedQuery = "";
	gaplessGermline = "";
	gaplessQuery = "";
	gapPenalties = 0;
	bool delOpen = false, insOpen = false;

	int k = germLen + seqLen - 1;
	int s = k, g = k;

	alignedGerm.resize(k + 1);
	alignedQuery.resize(k + 1);
	gaplessGermline.resize(k + 1);
	gaplessQuery.resize(k + 1);

	int seq_pos = final_posIJ[0];
	int germ_pos = final_posIJ[1];
	int next_matrix = DIAG, current_matrix, prev_matrix = DIAG;//always start at the diag position
	seqEnd = seq_pos - 1;
	germEnd = germ_pos - 1;
	numSeqIns[0] = 0;
	numSeqDel[0] = 0;
	numSeqDel[1] = 0;
	numSeqIns[1] = 0;

	int prevMatch = 0;
	int prevMisMatch = 0;

	double prevPenalty = 0;

	tracebackData[germ_pos][MATCHES] = 0;
	tracebackData[germ_pos][MISMATCHES] = 0;
	tracebackData[germ_pos][INDEL] = 0;
	tracebackData[germ_pos][GAPPENALTY] = 0;
	tracebackData[germ_pos][SCORE] = 0;


	endSeq = query.substr(seq_pos, query.length() - seq_pos);
	bool beg_of_seq = false;

	do{
		current_matrix = next_matrix;

		next_matrix = alignmentPath[current_matrix][seq_pos][germ_pos];

		if (next_matrix == START)
			break;

		switch (current_matrix){
		case DIAG:
			alignedGerm[k] = germline[germ_pos - 1];
			alignedQuery[k] = query[seq_pos - 1];
			gaplessGermline[g] = germline[germ_pos - 1];
			gaplessQuery[s] = query[seq_pos - 1];

			tracebackData[germ_pos - 1][SCORE] = scoreMatrix[next_matrix][seq_pos - 1][germ_pos - 1];
			tracebackData[germ_pos - 1][SEQPOS] = seq_pos - 1;
			tracebackData[germ_pos - 1][INDEL] = totalIndel;
			tracebackData[germ_pos - 1][ALGNPOS] = k;

			tracebackData[germ_pos - 1][GAPPENALTY] = gapPenalties;

			prevPenalty = 0;

			if (scoreMatrix[DIAG][seq_pos][germ_pos] < scoreMatrix[next_matrix][seq_pos - 1][germ_pos - 1]){
				prevMisMatch += 1;
			}
			else if (scoreMatrix[DIAG][seq_pos][germ_pos] > scoreMatrix[next_matrix][seq_pos - 1][germ_pos - 1]){
				prevMatch += 1;
			}

			tracebackData[germ_pos - 1][MATCHES] = prevMatch;
			tracebackData[germ_pos - 1][MISMATCHES] = prevMisMatch;

			s--;
			g--;
			seq_pos--;
			germ_pos--;

			break;

		case RIGHT:
			alignedGerm[k] = germline[germ_pos - 1];
			alignedQuery[k] = '-';
			gaplessQuery[s] = NOT_NT; //add an X in the seuence, but do not add any letter ot aligned germline

			numSeqIns[0] += 1;


			if (prev_matrix == RIGHT){
				seqIns[numSeqIns[1] - 1][1] += 1; //add 1 to the "seq insertions" list
			}
			else{
				seqIns[numSeqIns[1]][1] = 1; //add 1 to the "seq insertions" list
				numSeqIns[1]++;
			}

			if (next_matrix != RIGHT){
				gapPenalties += swParams.swGapOpen;
				tracebackData[germ_pos - 1][SCORE] = scoreMatrix[RIGHT][seq_pos][germ_pos] + swParams.swGapOpen;
			}
			else{
				gapPenalties += swParams.swExtendGap;
				tracebackData[germ_pos - 1][SCORE] = scoreMatrix[RIGHT][seq_pos][germ_pos] + swParams.swExtendGap;
			}

			tracebackData[germ_pos - 1][ALGNPOS] = k;

			totalIndel += 1;// prevIndel;
			prevPenalty = 0;


			//tracebackData[germ_pos - 1][SCORE] = scoreMatrix[next_matrix][seq_pos-1][germ_pos];
			tracebackData[germ_pos - 1][SEQPOS] = seq_pos;
			tracebackData[germ_pos - 1][INDEL] = totalIndel;
			tracebackData[germ_pos - 1][GAPPENALTY] = gapPenalties;
			tracebackData[germ_pos - 1][MATCHES] = prevMatch;
			tracebackData[germ_pos - 1][MISMATCHES] = prevMisMatch;


			seqIns[numSeqIns[1] - 1][0] = seq_pos - 1;
			germ_pos--;
			s--;
			break;

		case DOWN:
			alignedGerm[k] = '-';
			gaplessGermline[g] = NOT_NT;
			alignedQuery[k] = query[seq_pos - 1];

			numSeqDel[0] += 1;
			totalIndel += 1;

			if (prev_matrix == DOWN){
				gapPenalties += swParams.swExtendGap;
				prevPenalty += swParams.swExtendGap;
				seqDel[numSeqDel[1] - 1][1] += 1;
			}
			else{
				gapPenalties += swParams.swGapOpen;
				prevPenalty += swParams.swGapOpen;
				//prevIndel = totalIndel;
				seqDel[numSeqDel[1]][1] = 1;
				numSeqDel[1]++;
				tracebackData[germ_pos - 1][GAPPENALTY] = prevPenalty;
				//tracebackData[germ_pos - 1][INDEL] = totalIndel;
				//tracebackData[germ_pos - 1][SEQPOS] = seq_pos-1;
			}

			seqDel[numSeqDel[1] - 1][0] = seq_pos - 1;
			seq_pos--;
			g--;
			break;
		}
		prev_matrix = current_matrix;
		k--;
	} while (next_matrix != START);
	alignedGerm.erase(0, k + 1);
	alignedQuery.erase(0, k + 1);

	seqStart = seq_pos;
	germStart = germ_pos;

	if (alignedQuery[0] == '-'){
		numSeqIns[0] -= seqIns[numSeqIns[1] - 1][1];
		totalIndel -= seqIns[numSeqIns[1] - 1][1];
		double scoreAdjust = swParams.swExtendGap*(seqIns[numSeqIns[1] - 1][1] - 1) + swParams.swGapOpen;
		gapPenalties -= scoreAdjust;
		maxScore += scoreAdjust;
		germStart += seqIns[numSeqIns[1] - 1][1];
		s += seqIns[numSeqIns[1] - 1][1];
		alignedQuery.erase(0, seqIns[numSeqIns[1] - 1][1]);
		k += seqIns[numSeqIns[1] - 1][1];
		alignedGerm.erase(0, seqIns[numSeqIns[1] - 1][1]);
		seqIns[numSeqIns[1] - 1][1] = 0;
		numSeqIns[1] -= 1;
	}

	if (alignedGerm[0] == '-'){
		numSeqDel[0] -= seqDel[numSeqDel[1] - 1][1];
		totalIndel -= seqIns[numSeqDel[1] - 1][1];
		double scoreAdjust = swParams.swExtendGap*(seqDel[numSeqDel[1] - 1][1] - 1) + swParams.swGapOpen;
		gapPenalties -= scoreAdjust;
		maxScore += scoreAdjust;
		seqStart += seqDel[numSeqDel[1] - 1][1];
		g += seqDel[numSeqDel[1] - 1][1];
		k += seqDel[numSeqDel[1] - 1][1];
		alignedQuery.erase(0, seqDel[numSeqDel[1] - 1][1]);
		alignedGerm.erase(0, seqDel[numSeqDel[1] - 1][1]);
		seqDel[numSeqDel[1] - 1][1] = 0;
		numSeqDel[1] -= 1;
	}

	tracebackData[germStart][SCORE] = 0;
	tracebackData[germEnd + 1][SCORE] = maxScore;
	totalMatch = tracebackData[germStart][MATCHES];
	totalMismatch = tracebackData[germStart][MISMATCHES];

	for (int i = germStart; i <= germEnd + 1; i++){ //reverse the total indels/indel info at each index (remember we worked backwards)		
		tracebackData[i][MATCHES] = -1 * (tracebackData[i][MATCHES] - totalMatch);// +(tracebackData[i][MATCHES] - tracebackData[i + 1][MATCHES]);
		tracebackData[i][MISMATCHES] = -1 * (tracebackData[i][MISMATCHES] - totalMismatch);// +(tracebackData[i][MISMATCHES] - tracebackData[i + 1][MISMATCHES]);
		tracebackData[i][GAPPENALTY] = -1 * (tracebackData[i][GAPPENALTY] - gapPenalties);// +(tracebackData[i][GAPPENALTY] - tracebackData[i + 1][GAPPENALTY]);
		tracebackData[i][INDEL] = -1 * (tracebackData[i][INDEL] - totalIndel);// +(tracebackData[i][INDEL] - tracebackData[i + 1][INDEL]);
		tracebackData[i][ALGNPOS] -= k + 1;
	}

	begSeq = query.substr(0, seqStart);

	gaplessGermline.erase(0, g + 1);
	gaplessQuery.erase(0, s + 1);
	gaplessQuery = begSeq + gaplessQuery + endSeq;
}

void SWAlignment::printAlignment(){
	cout << "The alignment starts at positions: " << germStart << ", and " << seqStart << "(germ start/query start)" << endl;
	cout << "The alignment ends at positions: " << germEnd << ", and " << seqEnd << "(germ end/query end)" << endl;
	
	int a = 0;
	for (int i = 0; i < alignedGerm.length(); i++){
		cout << alignedGerm[i] << setw(2.5);
	}
	cout << ":----->Germline" << endl;
	for (int i = 0; i < alignedGerm.length(); i++){
		cout << alignedQuery[i] << setw(2.5);
	}
	cout << ":----->Query" << endl;
	for (int i = 0; i < alignedGerm.length(); i++){
		if (alignedGerm[i] != '-'){
			cout << germStart + a << setw(2.5);
			a += 1;
		}
		else{
			cout << "" << setw(2.5);
		}
	}
	cout << endl;

}
void SWAlignment::printTracebackAnnotation(int pos){
	cout << "RESULTS FOR POSITION " << pos << endl;
	cout << "Matches: " << tracebackData[pos][MATCHES] << endl;
	cout << "Mismatches: " << tracebackData[pos][MISMATCHES] << endl;
	cout << "GAPPENALTIES: " << tracebackData[pos][GAPPENALTY] << endl;
	cout << "INDELS: " << tracebackData[pos][INDEL] << endl;
	cout << "SCORE: " << tracebackData[pos][SCORE] << endl;
	cout << "SEQ POS: " << tracebackData[pos][SEQPOS] << endl;
	cout << "ALGN POS: " << tracebackData[pos][ALGNPOS] << endl;
	cout << "Results complete" << endl;
}