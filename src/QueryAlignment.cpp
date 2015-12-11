#include "QueryAlignment.h"
using namespace std;

/*CONSTRUCTORS*/
QueryAlignment::QueryAlignment(string & query_file, const string & inputFileFormat, const structvars::AlignmentProgramSettings & settings){
	/***READ QUERY FILE, DETERMINE NUMBER OF SEQUENCES AND MAX SEQUENCE LENGTH****/
	//if the input file is a fasta file, then quickly convert it to a tab, also input file name will be changed via reference
	
	if (inputFileFormat == "FASTA"){
		printf("Calcuating number of sequences\n");
		queryInfo = readfiles::convertFASTAtoTAB(query_file);
	}
	else{
		printf("Calcuating number of sequences\n");
		queryInfo = readfiles::readSeqLengths(query_file.c_str(), "TAB");
		
	}

	algnSettings = settings;

	full_seq_alignment = NULL;
	list_of_fft_plans = NULL;
	peptide_seq_alignment = NULL;

	germlineHits = NULL;
	aligned_to_clusters = NULL;
	finalized_sw_alignments = NULL;
	germlineSequences = NULL;
	seqPosStore = NULL;
	seqPeptideList = NULL;

	results_with_gaps = NULL;
	results_no_gaps = NULL;
	best_cluster_so_far = NULL;

	//inverseRatio = 1 / scoreFoldRatio;
	maxGermlineLength = 0;
	maxQueryLength = 0;
	queryStart = 0;
	queryEnd = 0;
	germStart = 0;
	germEnd = 0;
	numClusters = 0;
	numGermlines = 0;
	numClusterHits = 0;
	bestSeqDir[0] = '+';
	bestSeqDir[1] = '-';
	numUniqueChainHits = 0;
	numUniqueLocusHits = 0;
	
	annotationFields = { "FR1_FLANK", "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "CDR3_FLANK" };	
	std::string outfilename = "debuggingalignments.txt";
	//debugfile.open(outfilename.c_str());
	t = 0;
	tO = 0;
}

QueryAlignment::QueryAlignment(const structvars::Fileinfo & query_file_info, const structvars::AlignmentProgramSettings & settings){
	queryInfo = query_file_info;
	algnSettings = settings;

	full_seq_alignment = NULL;
	list_of_fft_plans = NULL;
	peptide_seq_alignment = NULL;

	germlineHits = NULL;
	aligned_to_clusters = NULL;
	finalized_sw_alignments = NULL;
	germlineSequences = NULL;
	seqPosStore = NULL;
	seqPeptideList = NULL;

	results_with_gaps = NULL;
	results_no_gaps = NULL;
	best_cluster_so_far = NULL;
	bestSeqDir[0] = '+';
	bestSeqDir[1] = '-';
	//scoreFoldRatio = 2;
	//inverseRatio = 1 / scoreFoldRatio;
	maxGermlineLength = 0;
	maxQueryLength = 0;
	queryStart = 0;
	queryEnd = 0;
	germStart = 0;
	germEnd = 0;
	numClusters = 0;
	numGermlines = 0;
	numClusterHits = 0;
	t = 0;
	numUniqueChainHits = 0;
	numUniqueLocusHits = 0;
	tO = 0;
	std::string outfilename = "debuggingalignments.txt";
	annotationFields = { "FR1_FLANK", "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "CDR3_FLANK" };
	//debugfile.open(outfilename.c_str());
}

QueryAlignment::QueryAlignment(string & query_file, const string & inputFileFormat,const string & germline_db_file, const string & cluster_file, const structvars::AlignmentProgramSettings & settings){
	/***READ QUERY FILE, DETERMINE NUMBER OF SEQUENCES AND MAX SEQUENCE LENGTH****/
	//if the input file is a fasta file, then quickly convert it to a tab, also input file name will be changed via reference
	if (inputFileFormat == "FASTA"){
		printf("Calcuating number of sequences\n");
		queryInfo = readfiles::convertFASTAtoTAB(query_file);
	}
	else{
		printf("Calcuating number of sequences\n");
		queryInfo = readfiles::readSeqLengths(query_file.c_str(), "TAB");
	}
	
	algnSettings = settings;
	full_seq_alignment = NULL;
	list_of_fft_plans = NULL;
	peptide_seq_alignment = NULL;

	germlineHits = NULL;
	aligned_to_clusters = NULL;
	finalized_sw_alignments = NULL;
	germlineSequences = NULL;
	seqPosStore = NULL;
	seqPeptideList = NULL;

	results_with_gaps = NULL;
	results_no_gaps = NULL;
	best_cluster_so_far = NULL;
	bestSeqDir[0] = '+';
	bestSeqDir[1] = '-';
	//scoreFoldRatio = 2;
	//inverseRatio = 1 / scoreFoldRatio;
	maxGermlineLength = 0;
	maxQueryLength = 0;
	queryStart = 0;
	queryEnd = 0;
	germStart = 0;
	germEnd = 0;
	numClusters = 0;
	numGermlines = 0;
	numUniqueChainHits = 0;
	numUniqueLocusHits = 0;
	numClusterHits = 0;
	ReadGermlineDatabase(germline_db_file, cluster_file);
	t = 0;
	tO = 0;
	std::string outfilename = "debuggingalignments.txt";
	annotationFields = { "FR1_FLANK", "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "CDR3_FLANK" };
	// debugfile.open(outfilename.c_str());
}

QueryAlignment::QueryAlignment(const structvars::Fileinfo & query_file_info, const string & germline_db_file, const string & cluster_file, const structvars::AlignmentProgramSettings & settings){
	/***READ QUERY FILE, DETERMINE NUMBER OF SEQUENCES AND MAX SEQUENCE LENGTH****/
	//if the input file is a fasta file, then quickly convert it to a tab, also input file name will be changed via reference
	queryInfo = query_file_info;
	algnSettings = settings;
	full_seq_alignment = NULL;
	list_of_fft_plans = NULL;
	peptide_seq_alignment = NULL;

	germlineHits = NULL;
	aligned_to_clusters = NULL;
	finalized_sw_alignments = NULL;
	germlineSequences = NULL;
	seqPosStore = NULL;
	seqPeptideList = NULL;
	bestSeqDir[0] = '+';
	bestSeqDir[1] = '-';
	results_with_gaps = NULL;
	results_no_gaps = NULL;
	best_cluster_so_far = NULL;

	//scoreFoldRatio = 2;
	//inverseRatio = 1 / scoreFoldRatio;
	maxGermlineLength = 0;
	maxQueryLength = 0;
	queryStart = 0;
	queryEnd = 0;
	germStart = 0;
	germEnd = 0;
	numClusters = 0;
	numGermlines = 0;
	numClusterHits = 0;
	numUniqueChainHits = 0;
	numUniqueLocusHits = 0;
	ReadGermlineDatabase(germline_db_file, cluster_file);
	t = 0;
	tO = 0;
	std::string outfilename = "debuggingalignments.txt";
	annotationFields = { "FR1_FLANK", "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "CDR3_FLANK" };
	//debugfile.open(outfilename.c_str());
}

QueryAlignment::~QueryAlignment(){

	if (list_of_fft_plans != NULL){
		full_seq_alignment = NULL;
		full_seq_alignment_RC = NULL;

		for (int i = 0; i < numFFTPlans; i++){
			for (int j = 0; j < 2; j++)
				delete list_of_fft_plans[i][j];
			delete[] list_of_fft_plans[i];			
		}
		delete[] list_of_fft_plans;		
	}
	
	for (int i = 0; i < numClusters; i++){
		aligned_to_clusters[i].queryInfoPointer = NULL;
	}

	if (peptide_seq_alignment != NULL)
		delete peptide_seq_alignment;

	if (results_no_gaps != NULL)
		delete[] results_no_gaps;
	if (results_with_gaps != NULL)
		delete[] results_with_gaps;

	DeleteGermlineVars();
	DeleteClusterVars();
	DeleteSeqPeptideVars();
	//debugfile.close();
}
/*END of CONSTRUCTORS*/

/*PUBLIC FUNCTIONS*/

/*ACCESSORS*/
structvars::Fileinfo QueryAlignment::ReturnQueryFileInfo(){
	return queryInfo;
}
structvars::Fileinfo QueryAlignment::GermlineQueryFileInfo(){
	return germlineSeqInfo;
}
structvars::Fileinfo QueryAlignment::ClusterQueryFileInfo(){
	return clusterSeqInfo;
}

/*functions*/
void QueryAlignment::ReadGermlineDatabase(const string & germline_file_name, const string & cluster_file_name,bool allowSubstrings){
	//NO LONGER A FUNCTIONAL FUNCTION//

	DeleteGermlineVars(); //if its not the frist time we called this function, then we need to delete the germline database we have made previously

	map<string, int> unique_locus_names;
	map<string, int> unique_chain_names;

	ReadGermlineFile(germline_file_name, unique_locus_names, unique_chain_names); //first read in all of the germlines defined in the germline file name

	//Now we know the size of the database, so we can start initlalizing some of the germline variables
	germlineHits = new structvars::GermlineAlignmentResults[numGermlines];
	//Initialize the size of the variables for storing germline results
	aligned_to_germline_scores.resize(numGermlines);

	germlineSequences = new string[numGermlines];
	for (int i = 0; i < numGermlines; i++){
		germlineSequences[i] = germlineDatabase[i].germlineSeqData.seq;
		aligned_to_germline_scores[i].resize(2);
	}
	
	//next read in the information where germline sequences have been pre-grouped into clusters
	printf("\tLoaded %lu germline sequences...\n", germlineDatabase.size());
	printf("\tReading in clusters from provided file...\n");
	ReadInClusterDatabase(cluster_file_name);
	printf("\tMapping pairwise cluster alignments...\n");
	AlignPairwiseClusters();	
	FindGermlineStart(); 
	printf("\tGermline sequences have been grouped into clusters.\n");
	printf("\tThere are a total number of %i clusters\n", numClusters);


	InitializeTransformVariables(allowSubstrings);
}

void QueryAlignment::ReadGermlineDatabase(const std::vector<std::string> & germline_files, int method, bool allowSubstrings){
	// DeleteGermlineVars(); //if its not the frist time we called this function, then we need to delete the germline database we have made previously

	map<string, int> unique_locus_names;
	map<string, int> unique_chain_names;

	for (int i = 0; i < germline_files.size(); i++){
		printf("%s\n", germline_files[i].c_str());
		ReadGermlineFile(germline_files[i], unique_locus_names, unique_chain_names); //first read in all of the germlines defined in the germline file name
	}

	//Now we know the size of the database, so we can start initlalizing some of the germline variables
	germlineHits = new structvars::GermlineAlignmentResults[numGermlines];
	aligned_to_germline_scores.resize(numGermlines);
	germlineSequences = new string[numGermlines];
	for (int i = 0; i < numGermlines; i++){
		germlineSequences[i] = germlineDatabase[i].germlineSeqData.seq;
		aligned_to_germline_scores[i].resize(2);
	}
	printf("\tLoaded %lu germline sequences...\n", germlineDatabase.size());
	printf("\tGrouping highly similar germline sequences into clusters...\n");
	//next read in the information concerning the clusters
	if (method == COPY)
		CopyGermlineDBToClusterDB(); //dont cluster, just make each germline an indivdual cluster result 
	else if (method == MAKE)
		GroupGermlinesIntoClusters();  

	printf("\tMapping pairwise cluster alignments...\n");
	AlignPairwiseClusters();
	FindGermlineStart();
	printf("\tGermline sequences have been grouped into clusters.\n");
	printf("\tThere are a total number of %i clusters\n", numClusters);
	InitializeTransformVariables(allowSubstrings);
}

//this function will take the germline databse and group it into clusters based on a predefined cutoff
void QueryAlignment::GroupGermlinesIntoClusters(){
	vector<vector<int>> clusterDB;
	//vector<int> tempV;
	
	int alignLen = 2*germlineSeqInfo.maxSeqLen;
	vector<int> germlineModel;	
	vector<int> clusteredIndex(numGermlines);
	germlineModel.resize(alignLen);
	for (int i = 0; i < germlineSeqInfo.maxSeqLen; i++)
		germlineModel[i] = i + 1;
	for (int i = germlineSeqInfo.maxSeqLen; i < alignLen; i++)
		germlineModel[i] = 0;
	
	FFTAlign gapplessAlignment(alignLen, germlineModel,germlineModel, algnSettings.fftParams);


	complex_t **germlinesComplex, **germlinesFFT, *temp;
	germlinesComplex = (complex_t**)malloc(numGermlines*sizeof(complex_t*));	
	germlinesFFT = (complex_t**)malloc(numGermlines*sizeof(complex_t*));
	temp = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
	double score,perID;
	int diag;
	int avgAlgnLen;
	
	int numAlgnRows = numGermlines*(numGermlines - 1) / 2, si, sj, ei , ej;
	
	vector<vector<double>> pairwiseAlignments(numAlgnRows, { vector<double>(8) });

	/*convert all germline sequences to FFT*/
	for (int i = 0; i < numGermlines; i++){		
		germlinesComplex[i] = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
		germlinesFFT[i] = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));		
		ffthelper::convertSeq(temp, germlineDatabase[i].germlineSeqData.seq); //convert string to complex numbers	
		ffthelper::PadSeqZerosEnd(temp, germlinesComplex[i], germlineDatabase[i].germlineSeqData.seqLength, gapplessAlignment.alignmentLength);
		gapplessAlignment.ForwardFFT(germlinesComplex[i], germlinesFFT[i]);
	}

	int row = 0;
	int pairAlgnLen;
	for (int i = 0; i < numGermlines; i++){
		for (int j = i + 1; j < numGermlines; j++){
			if (germlineDatabase[i].locus != germlineDatabase[j].locus){//ensure that we only cluster together germlines from the same locus
				pairwiseAlignments[row][0] = i;
				pairwiseAlignments[row][1] = j;
				pairwiseAlignments[row][2] = -1;
				pairwiseAlignments[row][3] = -1;				
				pairwiseAlignments[row][4] = 0;
				pairwiseAlignments[row][5] = 0;
				pairwiseAlignments[row][6] = 0;
				pairwiseAlignments[row][7] = 0;
			}
			else{
				for (int k = 0; k < gapplessAlignment.alignmentLength; k++){
					temp[k][0] = germlinesFFT[j][k][0] / gapplessAlignment.alignmentLength;
					temp[k][1] = -1 * germlinesFFT[j][k][1] / gapplessAlignment.alignmentLength;
				}
				gapplessAlignment.AlignSequences(germlinesFFT[i], temp);
				score = gapplessAlignment.bestDiagScore;
				diag = gapplessAlignment.bestDiag;

				if (abs(1 - ((germlineDatabase[i].germlineSeqData.seqLength*1.000) / (1.000*germlineDatabase[j].germlineSeqData.seqLength))) <= 0.2)
					avgAlgnLen = germlineDatabase[i].germlineSeqData.seqLength + germlineDatabase[j].germlineSeqData.seqLength;
				else
					avgAlgnLen = 2 * min(germlineDatabase[i].germlineSeqData.seqLength, germlineDatabase[j].germlineSeqData.seqLength);

				perID = 2 * score / avgAlgnLen;
				pairwiseAlignments[row][0] = i;
				pairwiseAlignments[row][1] = j;
				pairwiseAlignments[row][2] = score;
				pairwiseAlignments[row][3] = perID;
				si = gapplessAlignment.alignmentModel[diag][1] - 1;//starting base of alignment for first sequence
				sj = gapplessAlignment.alignmentModel[diag][3] - 1;//starting base of alingment for second sequence
				pairAlgnLen = min(gapplessAlignment.alignmentModel[diag][5], germlineDatabase[i].germlineSeqData.seqLength - si);
				pairAlgnLen = min(pairAlgnLen, germlineDatabase[j].germlineSeqData.seqLength - sj);
				ei = si + pairAlgnLen - 1;// min(alignmentInfo[diag][2] - 1, germlineDatabase[i].germlineSeqData.seqLength - 1);
				ej = sj + pairAlgnLen - 1;// min(alignmentInfo[diag][4] - 1, germlineDatabase[j].germlineSeqData.seqLength - 1);
				pairwiseAlignments[row][4] = si;
				pairwiseAlignments[row][5] = sj;
				pairwiseAlignments[row][6] = ei;
				pairwiseAlignments[row][7] = ej;
			}
			
			row++;
		}
	}

	//sort(pairwiseAlignments.begin(), pairwiseAlignments.end(), [](const vector<double> & a, const vector<double> & b){ return (a[2] > b[2]); }); //sort 2d vector using the second element in the vector (takes advantage of lamda functions)
		
	SWAlignment pairAlgn(germlineSeqInfo.maxSeqLen, germlineSeqInfo.maxSeqLen, algnSettings.swParams);
	int unclusteredSeqs = 0;
	int g1, g2, tempI1, tempI2;
	bool no_more_clusters = false, no_gaps_in_other_members;
	int not_clustered, current_count, current_cluster;
	
	//ow = (tempI1*(numGermlines - 1) - (tempI1*(tempI1 - 1)) / 2) + (tempI2 - (tempI1 + 1)); //for mapping x and y alignments to row in parisie alignments
	unclusteredSeqs = 0;
	current_count = 0;
	not_clustered = 0;	
	int numSWALIGNED = 0;
	for (int k = 0; k<numAlgnRows; k++){
		g1 = int(pairwiseAlignments[k][0]);
		g2 = int(pairwiseAlignments[k][1]);
		current_count++;
		
		if (pairwiseAlignments[k][3] >= algnSettings.clusterGermlineCutoff && germlineDatabase[pairwiseAlignments[k][0]].locus == germlineDatabase[pairwiseAlignments[k][1]].locus){ //the pariwise alignmetn is above the expected alignment				
		
			if (clusteredIndex[g1] == 0 && clusteredIndex[g2] == 0){
				clusterDB.push_back(vector<int>{g1, g2});//form a new cluster
				clusteredIndex[g1] = clusterDB.size();
				clusteredIndex[g2] = clusterDB.size();
			}
			else if (clusteredIndex[g1] > 0 && clusteredIndex[g2] == 0){ //consider adding G2 to cluster containing G1				
				current_cluster = clusteredIndex[g1] - 1;
				clusterDB[current_cluster].push_back(g2);
				clusteredIndex[g2] = clusteredIndex[g1];				
			}
			else if (clusteredIndex[g2] > 0 && clusteredIndex[g1] == 0){ //consider adding G1 to cluster containing g2
				current_cluster = clusteredIndex[g2] - 1;				
				clusterDB[current_cluster].push_back(g1);
				clusteredIndex[g1] = clusteredIndex[g2];
				
			}
			else if (clusteredIndex[g2] != clusteredIndex[g1]){ //consider mergine together the clusters which contain g1 and g2
				tempI1 = clusteredIndex[g2] - 1;
				for (int l = 0; l < clusterDB[tempI1].size(); l++){
					tempI2 = clusterDB[tempI1][l];
					clusteredIndex[tempI2] = clusteredIndex[g1];
				}
				clusterDB[clusteredIndex[g1] - 1].insert(clusterDB[clusteredIndex[g1] - 1].end(), clusterDB[tempI1].begin(), clusterDB[tempI1].end());
				clusterDB[tempI1].clear();
				
			}
		}		
	}

	int i = 0;
	do{
		if (clusterDB[i].size() == 0)
			clusterDB.erase(clusterDB.begin() + i);
		else{
			sort(clusterDB[i].begin(), clusterDB[i].end(), [](int a, int b){ return (b>a); });
			i++;
		}
	} while (i < clusterDB.size());

	for (int i = 0; i<numGermlines; i++){ //remaining germlines not placed into clusters are singletons
		if (clusteredIndex[i] == 0){
			clusterDB.push_back(vector<int>{i});
			clusteredIndex[i] = clusterDB.size();
		}
	}
	vector<vector<int>> finalizedClusters;
	vector<string> finalizedConsSeq;
	for (int i = 0; i < clusterDB.size(); i++){
		if (clusterDB[i].size()>1)
			FindConsensusSeqInCluster(clusterDB[i], pairwiseAlignments, finalizedClusters, finalizedConsSeq);
		else{
			finalizedClusters.push_back(clusterDB[i]);
			finalizedConsSeq.push_back(germlineDatabase[clusterDB[i][0]].germlineSeqData.seq);
		}

	}
	
	/*
	unclusteredSeqs = 0;
	current_count = 0;
	not_clustered = 0;
	cout << numAlgnRows << endl;
	int numSWALIGNED = 0;
	for (int k = 0; k<numAlgnRows; k++){
		g1 = int(pairwiseAlignments[k][0]);
		g2 = int(pairwiseAlignments[k][1]);
		current_count++;
		if (pairwiseAlignments[k][3] >= algnSettings.clusterGermlineCutoff){ //the pariwise alignmetn is above the expected alignment				
			if (!aligned[g1][g2]){
				pairAlgn.OverlapAlignComplete(germlineDatabase[g1].germlineSeqData.seq, germlineDatabase[g2].germlineSeqData.seq); //align each cluster to the "seed cluster consensus". 
				hasGaps[g1][g2] = pairAlgn.numSeqIns[0] + pairAlgn.numSeqDel[0] > 0;
				hasGaps[g2][g1] = hasGaps[g1][g2];
				aligned[g1][g2] = true;
				aligned[g2][g1] = true;
				numSWALIGNED++;
			}
				
			if (!hasGaps[g1][g2]){
				if (clusteredIndex[g1] == 0 && clusteredIndex[g2] == 0){ 
					clusterDB.push_back(vector<int>{g1, g2});//form a new cluster
					clusteredIndex[g1] = clusterDB.size();
					clusteredIndex[g2] = clusterDB.size();
				}
				else if (clusteredIndex[g1] > 0 && clusteredIndex[g2] == 0){ //consider adding G2 to cluster containing G1
					no_gaps_in_other_members = true;
					current_cluster = clusteredIndex[g1] - 1;
					for (int j = 0; j < clusterDB[current_cluster].size(); j++){
						tempI1 = clusterDB[current_cluster][j];

						if (!aligned[tempI1][g2]){
							pairAlgn.OverlapAlignComplete(germlineDatabase[tempI1].germlineSeqData.seq, germlineDatabase[g2].germlineSeqData.seq); //align each cluster to the "seed cluster consensus". 
							hasGaps[tempI1][g2] = pairAlgn.numSeqIns[0] + pairAlgn.numSeqDel[0] > 0;
							hasGaps[g2][tempI1] = hasGaps[tempI1][g2];
							aligned[tempI1][g2] = true;
							aligned[g2][tempI1] = true;
							numSWALIGNED++;
						}
						if (hasGaps[tempI1][g2])
							no_gaps_in_other_members = false;
						
					}
					if (no_gaps_in_other_members){//yes add g2 to g1 cluster
						clusterDB[current_cluster].push_back(g2);
						clusteredIndex[g2] = clusteredIndex[g1];
					}
	
				}
				else if (clusteredIndex[g2] > 0 && clusteredIndex[g1] == 0){ //consider adding G1 to cluster containing g2
					no_gaps_in_other_members = true;
					current_cluster = clusteredIndex[g2] - 1;
					for (int j = 0; j < clusterDB[current_cluster].size(); j++){
						tempI2 = clusterDB[current_cluster][j];
						if (!aligned[tempI2][g1]){
							pairAlgn.OverlapAlignComplete(germlineDatabase[tempI2].germlineSeqData.seq, germlineDatabase[g1].germlineSeqData.seq); //align each cluster to the "seed cluster consensus". 
							hasGaps[tempI2][g1] = pairAlgn.numSeqIns[0] + pairAlgn.numSeqDel[0] > 0;							
							hasGaps[g1][tempI2] = hasGaps[tempI2][g1];
							aligned[tempI2][g1] = true;
							aligned[g1][tempI2] = true;
							numSWALIGNED++;
						}
						if (hasGaps[tempI2][g1])
							no_gaps_in_other_members = false;
					}
					if (no_gaps_in_other_members){//yes add g1 to g2 cluster
						clusterDB[current_cluster].push_back(g1);
						clusteredIndex[g1] = clusteredIndex[g2];
					}
					
				}
				else if (clusteredIndex[g2] != clusteredIndex[g1]){ //consider mergine together the clusters which contain g1 and g2
					no_gaps_in_other_members = true;
					current_cluster = clusteredIndex[g2] - 1;
					for (int j = 0; j < clusterDB[current_cluster].size(); j++){
						tempI2 = clusterDB[current_cluster][j];
						if (!aligned[tempI2][g1]){
							pairAlgn.OverlapAlignComplete(germlineDatabase[tempI2].germlineSeqData.seq, germlineDatabase[g1].germlineSeqData.seq); //align each cluster to the "seed cluster consensus". 
							hasGaps[tempI2][g1] = pairAlgn.numSeqIns[0] + pairAlgn.numSeqDel[0] > 0;
							hasGaps[g1][tempI2] = hasGaps[tempI2][g1];
							aligned[tempI2][g1] = true;
							aligned[g1][tempI2] = true;
							numSWALIGNED++;
						}
						if (hasGaps[tempI2][g1])
							no_gaps_in_other_members = false;
					}
					current_cluster = clusteredIndex[g1] - 1;
					for (int j = 0; j < clusterDB[current_cluster].size(); j++){
						tempI1 = clusterDB[current_cluster][j];
						if (!aligned[tempI1][g2]){
							pairAlgn.OverlapAlignComplete(germlineDatabase[tempI1].germlineSeqData.seq, germlineDatabase[g2].germlineSeqData.seq); //align each cluster to the "seed cluster consensus". 
							hasGaps[tempI1][g2] = pairAlgn.numSeqIns[0] + pairAlgn.numSeqDel[0] > 0;
							hasGaps[g2][tempI1] = hasGaps[tempI1][g2];
							aligned[tempI1][g2] = true;
							aligned[g2][tempI1] = true;
							numSWALIGNED++;
						}
						if (hasGaps[tempI1][g2])
							no_gaps_in_other_members = false;
					}

					if (no_gaps_in_other_members){
						tempI1 = clusteredIndex[g2] - 1;
						for (int l = 0; l < clusterDB[tempI1].size(); l++){
							tempI2 = clusterDB[tempI1][l];
							clusteredIndex[tempI2] = clusteredIndex[g1];
						}
						clusterDB[clusteredIndex[g1] - 1].insert(clusterDB[clusteredIndex[g1] - 1].end(), clusterDB[tempI1].begin(), clusterDB[tempI1].end());
						clusterDB[tempI1].clear();
					}				
				}
			}					
		}
		cout << k << ": " << numSWALIGNED << ";  ";
	}

	int i = 0;
	do{
		if (clusterDB[i].size() == 0)
			clusterDB.erase(clusterDB.begin() + i);
		else{
			sort(clusterDB[i].begin(), clusterDB[i].end(), [](int a, int b){ return (b>a); }); 
			i++;
		}
	} while (i < clusterDB.size());

	for (int i = 0; i<numGermlines; i++){ //remaining germlines not placed into clusters are singletons
		if (clusteredIndex[i] == 0){
			clusterDB.push_back(vector<int>{i});
			clusteredIndex[i] = clusterDB.size();
		}
	}
	*/

	//OK we have now formed our clusters, lets go ahead and update our germline cluster info now that we know how many clusters to expect
	numClusters = finalizedClusters.size();
	aligned_to_clusters = new GermlineCluster[numClusters];
	int minSLen = 0, maxSLen = 0;
	int totalLen = 0;
	for (int i = 0; i < finalizedClusters.size(); i++){		
		aligned_to_clusters[i].consensus = finalizedConsSeq[i];
		totalLen += aligned_to_clusters[i].consensus.length();
		aligned_to_clusters[i].clusterSeqLength = finalizedConsSeq[i].length();
		if (aligned_to_clusters[i].clusterSeqLength > maxSLen)
			maxSLen = aligned_to_clusters[i].clusterSeqLength;
		if (aligned_to_clusters[i].clusterSeqLength<minSLen)
			minSLen = aligned_to_clusters[i].clusterSeqLength;
		aligned_to_clusters[i].clusterIndex = i;
		aligned_to_clusters[i].AddGermlineMemberInfo(finalizedClusters[i], germlineDatabase);		
		aligned_to_clusters[i].algnParams = algnSettings;
	}

	clusterSeqInfo.numSeqs = numClusters;
	clusterSeqInfo.maxSeqLen = maxSLen;
	clusterSeqInfo.minSeqLen = minSLen;
	clusterSeqInfo.avgSeqLen = (totalLen / numClusters) + 1;
	maxClusterLength = clusterSeqInfo.maxSeqLen; //"max germline length" just means the maximum length of the database we use to run FFT, so in fact its max cluster length

	for (int i = 0; i < numGermlines; i++){
		std::free(germlinesComplex[i]);
		std::free(germlinesFFT[i]);
	}

	std::free(germlinesComplex);
	std::free(germlinesFFT);
	std::free(temp);	
}

void QueryAlignment::FindConsensusSeqInCluster(vector<int> & cluster, const vector<vector<double>> & pairwiseAlignmentScores, vector<vector<int>> & finalClusterMembers, vector<string> & consSeqs){	
	
	vector<vector<int>> clusterData(cluster.size(), { vector<int>(2) });
	vector<string> gaplessClusterConsensus;
	vector<string> clusterSeeds;
	vector<vector<string>> clusterSeedMemberSeqs;
	vector<vector<int>> clusterMembers;	
		
	int maxLen = 0;
	for (int i = 0; i < cluster.size(); i++){
		clusterData[i][0] = cluster[i];
		clusterData[i][1] =  germlineDatabase[cluster[i]].seqLength;
		if (germlineDatabase[cluster[i]].seqLength > maxLen)
			maxLen = germlineDatabase[cluster[i]].seqLength;
	}
		
	std::sort(clusterData.begin(), clusterData.end(), [](const vector<int> & a, const vector<int> & b){ return (a[1] > b[1]);}); //sort 2d vector using the second element in the vector (takes advantage of lamda functions)//sort clusters by descending length
	
	/*determine the maximum consensus sequence length assuming there are no gaps*/
	int x, y,row;
	x = clusterData[0][0];
	double minP=0, maxP=germlineDatabase[x].germlineSeqData.seqLength-1;
	for (int i = 1; i < clusterData.size(); i++){
		y = clusterData[i][0];
		if (x < y){
			row = (x*(numGermlines - 1) - (x*(x- 1)) / 2) + (y - (x + 1)); //for mapping x and y alignments to row in parisie alignments
			minP = min(pairwiseAlignmentScores[row][4] - pairwiseAlignmentScores[row][5], minP);
			maxP = max(pairwiseAlignmentScores[row][4] + (germlineDatabase[y].germlineSeqData.seqLength-1-pairwiseAlignmentScores[row][5]), maxP);			
		}
		else{
			row = (y*(numGermlines - 1) - (y*(y - 1)) / 2) + (x - (y + 1)); //for mapping x and y alignments to row in parisie alignments
			minP = min(pairwiseAlignmentScores[row][5] - pairwiseAlignmentScores[row][4], minP);
			maxP = max(pairwiseAlignmentScores[row][5] + (germlineDatabase[y].germlineSeqData.seqLength - 1 - pairwiseAlignmentScores[row][4]), maxP);
		}		
	}
	/**/

	int front = abs(minP),to_insert,to_append;
	int end = maxP - (germlineDatabase[x].germlineSeqData.seqLength - 1);

	gaplessClusterConsensus.push_back(germlineDatabase[x].germlineSeqData.seq);
	gaplessClusterConsensus[0].insert(0, front, NOT_NT);	
	gaplessClusterConsensus[0].append(end, NOT_NT);
	
	/*adjust all sequences to the same length*/
	for (int i = 1; i < clusterData.size(); i++){		
		y = clusterData[i][0];
		gaplessClusterConsensus.push_back(germlineDatabase[y].germlineSeqData.seq);
		if (x < y){
			row = (x*(numGermlines - 1) - (x*(x - 1)) / 2) + (y - (x + 1)); //for mapping x and y alignments to row in parisie alignments
			to_insert = pairwiseAlignmentScores[row][4] - pairwiseAlignmentScores[row][5] + front;			
			gaplessClusterConsensus[i].insert(0,to_insert, NOT_NT);	
			to_append = gaplessClusterConsensus[0].length() - gaplessClusterConsensus[i].length();
			gaplessClusterConsensus[i].append(to_append, NOT_NT);
		}
		else{
			row = (y*(numGermlines - 1) - (y*(y - 1)) / 2) + (x - (y + 1)); //for mapping x and y alignments to row in parisie alignments
			to_insert = pairwiseAlignmentScores[row][5] - pairwiseAlignmentScores[row][4] + front;
			gaplessClusterConsensus[i].insert(0, to_insert, NOT_NT);
			to_append = gaplessClusterConsensus[0].length() - gaplessClusterConsensus[i].length();
			gaplessClusterConsensus[i].append(to_append, NOT_NT);			
		}
	}
	/****/
	
	string consSeq;
	consSeq.resize(gaplessClusterConsensus[0].length());
	string options = "ACGT";
	vector<vector<int>> basePrevalence(4, { vector<int>(2) });
	
	for (int i = 0; i < consSeq.length(); i++){		
		for (int j = 0; j < 4; j++){
			basePrevalence[j][0] = j;
			basePrevalence[j][1] = 0;
		}
		for (int j = 0; j < gaplessClusterConsensus.size(); j++){ //cacluate base consensus at position
			switch (gaplessClusterConsensus[j][i]){
			case 'A':
				basePrevalence[0][1]++;
				break;
			case 'C':
				basePrevalence[1][1]++;
				break;
			case 'G':
				basePrevalence[2][1]++;
				break;
			case 'T':
				basePrevalence[3][1]++;
				break;
			}
		}
		sort(basePrevalence.begin(), basePrevalence.end(), [](const vector<int> & a, const vector<int> & b){ return (a[1] > b[1]); }); //sort 2d vector using the second element in the vector (takes advantage of lamda functions)//sort clusters by descending length
		consSeq[i] = basePrevalence[0][1] > 0 ? options[basePrevalence[0][0]] : NOT_NT;
		
	}

	/*now we have a potentail cluster seed. use this consensus to align all sequences to the cluster seed and determine if a gap is found*/
	//clusterSeeds.push_back(vector < string > {consSeq});
	clusterSeeds.push_back(consSeq);
	clusterSeedMemberSeqs.resize(1);
	SWAlignment pairAlgn(consSeq.length(), consSeq.length(), algnSettings.swParams);
	bool hasGaps, allClustersHaveGaps;
	int clus_with_gaps;
	for (int i = 0; i < clusterData.size(); i++){
		hasGaps = false;
		clus_with_gaps = 0;
		for (int j = 0; j < clusterSeeds.size(); j++){
			pairAlgn.OverlapAlignComplete(clusterSeeds[j], gaplessClusterConsensus[i]); //align each cluster to the "seed cluster consensus". 
			if (pairAlgn.numSeqIns[0] + pairAlgn.numSeqDel[0] > 0)
				clus_with_gaps++;
			else{
				clusterSeedMemberSeqs[j].push_back(gaplessClusterConsensus[i]);
			}
		}
		if (clus_with_gaps == clusterSeeds.size()){ //gaps were found in all sequences
			clusterSeeds.push_back(gaplessClusterConsensus[i]); //add to the cluster consensus sequences
			clusterSeedMemberSeqs.push_back(vector < string > {gaplessClusterConsensus[i]}); //add to the list of sequence members in that new cluster			
		}
	}

	clusterMembers.resize(clusterSeeds.size());
	if (clusterSeeds.size() > 1){
		/*we now have a set of seeded clusters such that every member within the starting cluster should align to at least one cluster with no gaps*/
		for (int k = 0; k < clusterSeeds.size(); k++){
			for (int i = 0; i < clusterSeeds[k].length(); i++){
				for (int j = 0; j < 4; j++){
					basePrevalence[j][0] = j;
					basePrevalence[j][1] = 0;
				}

				for (int j = 0; j < clusterSeedMemberSeqs[k].size(); j++){ //cacluate base consensus at position
					switch (clusterSeedMemberSeqs[k][j][i]){
					case 'A':
						basePrevalence[0][1]++;
						break;
					case 'C':
						basePrevalence[1][1]++;
						break;
					case 'G':
						basePrevalence[2][1]++;
						break;
					case 'T':
						basePrevalence[3][1]++;
						break;
					}
				}
				std::sort(basePrevalence.begin(), basePrevalence.end(), [](const vector<int> & a, const vector<int> & b){ return (a[1] > b[1]); }); //sort 2d vector using the second element in the vector (takes advantage of lamda functions)//sort clusters by descending length
				clusterSeeds[k][i] = basePrevalence[0][1] > 0 ? options[basePrevalence[0][0]] : NOT_NT;			 
			}
		}

		/*now we want to make sure each member is aligned to the best possible cluster seed*/
		double bestScore;
		int bestCluster;
		for (int i = 0; i < clusterSeedMemberSeqs.size(); i++){
			clusterSeedMemberSeqs[i].clear();
		}
		for (int i = 0; i < clusterData.size(); i++){
			bestScore = std::numeric_limits<double>::lowest();;
			bestCluster = -1;
			for (int j = 0; j < clusterSeeds.size(); j++){
				pairAlgn.OverlapAlignComplete(clusterSeeds[j], gaplessClusterConsensus[i]); //align each cluster to the "seed cluster consensus". 
				if (pairAlgn.numSeqIns[0] + pairAlgn.numSeqDel[0] == 0 && pairAlgn.maxScore > bestScore){
					bestScore = pairAlgn.maxScore;
					bestCluster = j;
				}
			}
			if (bestCluster > -1){
				clusterSeedMemberSeqs[bestCluster].push_back(gaplessClusterConsensus[i]);
				clusterMembers[bestCluster].push_back(clusterData[i][0]);
			}
			else{
				clusterSeeds.push_back(gaplessClusterConsensus[i]); //add a new consensus cluster sequence
				clusterSeedMemberSeqs.push_back(vector < string > {gaplessClusterConsensus[i]}); //add to the list of sequence members in that new cluster							
				clusterMembers.push_back(vector < int > {clusterData[i][0]});
			}
		}

		for (int k = 0; k < clusterSeeds.size(); k++){
			for (int i = 0; i < clusterSeeds[k].length(); i++){
				for (int j = 0; j < 4; j++){
					basePrevalence[j][0] = j;
					basePrevalence[j][1] = 0;
				}

				for (int j = 0; j < clusterSeedMemberSeqs[k].size(); j++){ //cacluate base consensus at position
					switch (clusterSeedMemberSeqs[k][j][i]){
					case 'A':
						basePrevalence[0][1]++;
						break;
					case 'C':
						basePrevalence[1][1]++;
						break;
					case 'G':
						basePrevalence[2][1]++;
						break;
					case 'T':
						basePrevalence[3][1]++;
						break;
					}
				}
				sort(basePrevalence.begin(), basePrevalence.end(), [](const vector<int> & a, const vector<int> & b){ return (a[1] > b[1]); }); //sort 2d vector using the second element in the vector (takes advantage of lamda functions)//sort clusters by descending length				
				clusterSeeds[k][i] = basePrevalence[0][1] > 0 ? options[basePrevalence[0][0]] : NOT_NT;
			}
		}
	}
	else{
		for (int i = 0; i < clusterData.size(); i++)
			clusterMembers[0].push_back(clusterData[i][0]);
	}
		
	for (int i = 0; i < clusterSeeds.size(); i++){
		clusterSeeds[i].erase(clusterSeeds[i].find_last_not_of(NOT_NT) + 1); //trim right
		clusterSeeds[i].erase(0,clusterSeeds[i].find_first_not_of(NOT_NT));//trim left 		
		finalClusterMembers.push_back(clusterMembers[i]);
		consSeqs.push_back(clusterSeeds[i]);
	}	
	//by this point, all new clusters should be made and there should be no gaps between members within a cluster
}

void QueryAlignment::AlignPairwiseClusters(){

	int alignLen = 2 * maxClusterLength;
	vector<int> germlineModel, germlineQueryModel;
	germlineModel.resize(alignLen);
	germlineQueryModel.resize(alignLen);

	for (int i = 0; i < maxClusterLength; i++){
		germlineModel[i] = i + 1;
		germlineQueryModel[i] = 0;
	}
	for (int i = maxClusterLength; i < alignLen; i++){
		germlineModel[i] = 0;
		germlineQueryModel[i] = i-maxClusterLength + 1;
	}
	
	FFTAlign gapplessAlignment(alignLen, germlineQueryModel, germlineModel, algnSettings.fftParams);
	SWAlignment pairAlgn(maxClusterLength, maxClusterLength, algnSettings.swParams);
	complex_t **clustersComplex,**clustersComplexQuery, **clustersFFT,**clustersFFTQuery, *temp;

	clustersComplex = (complex_t**)malloc(numClusters*sizeof(complex_t*));
	clustersComplexQuery = (complex_t**)malloc(numClusters*sizeof(complex_t*));
	clustersFFT = (complex_t**)malloc(numClusters*sizeof(complex_t*));
	clustersFFTQuery = (complex_t**)malloc(numClusters*sizeof(complex_t*));
	temp = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
	
	/*convert all clusters sequences to FFT*/
	for (int i = 0; i < numClusters; i++){
		clustersComplex[i] = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
		clustersComplexQuery[i] = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
		clustersFFT[i] = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
		clustersFFTQuery[i] = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
		ffthelper::convertSeq(temp, aligned_to_clusters[i].consensus); //convert string to complex numbers	
		ffthelper::PadSeqZerosEnd(temp, clustersComplex[i], aligned_to_clusters[i].clusterSeqLength, gapplessAlignment.alignmentLength);
		ffthelper::PadSeqZerosFront(temp, clustersComplexQuery[i], aligned_to_clusters[i].clusterSeqLength, gapplessAlignment.alignmentLength);		
		gapplessAlignment.ForwardFFT(clustersComplex[i], clustersFFT[i]);
		gapplessAlignment.ForwardFFT(clustersComplexQuery[i], clustersFFTQuery[i]);
		for (int k = 0; k < gapplessAlignment.alignmentLength; k++){
			clustersFFT[i][k][0] = clustersFFT[i][k][0] / gapplessAlignment.alignmentLength;
			clustersFFT[i][k][1] = -1*clustersFFT[i][k][1] / gapplessAlignment.alignmentLength;
		}
	}

	vector<vector<structvars::PairwiseClusterAlignment>> clusterPairAlgnInfo(numClusters, { vector<structvars::PairwiseClusterAlignment>(numClusters) });
	vector<vector<int>> oppositeInsAlgnInfo,oppositeDelAlgnInfo, insAlgnInfo, delAlgnInfo;
	int qPos, gPos, minD, maxD, currentD, shift, ins, del,numMut;
	//int coordinates[5];
	for (int i = 0; i < numClusters; i++){
		for (int j = i; j < numClusters; j++){						
			currentD = 0;
			ins = -2;
			del = -2;
			oppositeInsAlgnInfo.clear(); //these will store ins/del events of the oppostie alignment (i.e. align the "i" sequence to "j" sequence as compared to "j" to "i")
			oppositeDelAlgnInfo.clear(); //these will store ins/del events of the oppostie alignment (i.e. align the "i" sequence to "j" sequence as compared to "j" to "i")
			insAlgnInfo.clear();
			delAlgnInfo.clear();
			if (i != j){
				gapplessAlignment.AlignSequences(clustersFFTQuery[j], clustersFFT[i]); //find the best alignment diagonal
				gapplessAlignment.IdentifyAdditionalAlignments(aligned_to_clusters[i].clusterSeqLength); //search for additional fft hits (cause by gaps)
				clusterPairAlgnInfo[i][j].fftMutations = gapplessAlignment.numGaps;				
				clusterPairAlgnInfo[j][i].fftMutations = gapplessAlignment.numGaps;
				qPos = (gapplessAlignment.alignmentModel[gapplessAlignment.bestDiag][1] - 1) - (maxClusterLength - aligned_to_clusters[j].clusterSeqLength);
				gPos = gapplessAlignment.alignmentModel[gapplessAlignment.bestDiag][3] - 1;
				clusterPairAlgnInfo[i][j].bestDiag = qPos - gPos;
				clusterPairAlgnInfo[j][i].bestDiag = gPos - qPos;
				
				pairAlgn.OverlapAlignComplete(aligned_to_clusters[j].consensus, aligned_to_clusters[i].consensus);
				qPos = pairAlgn.seqStart-1;
				gPos = pairAlgn.germStart-1;
				minD = clusterPairAlgnInfo[i][j].bestDiag;
				maxD = clusterPairAlgnInfo[i][j].bestDiag;
				
				for (int k = 0; k < pairAlgn.alignedQuery.length(); k++){
					if (pairAlgn.alignedGerm[k] != '-')
						gPos++;
					else{
						numMut = qPos - ins;
						if (numMut == 1){
							oppositeInsAlgnInfo[oppositeInsAlgnInfo.size() - 1][1]++;//insertion information for aligning the sequences in the REVERSE ORDER (i.e. aligning the i'th sequence to the j'th sequence)
						}
						else{
							oppositeInsAlgnInfo.push_back(vector<int>{gPos, 1});//insertion information for aligning the sequences in the REVERSE ORDER (i.e. aligning the i'th sequence to the j'th sequence)						
						}
						ins = qPos;
						delAlgnInfo.push_back(vector<int>{qPos+1, 1});//deletion information for aligning the j'th sequence to the i'th sequence)
					}
					if (pairAlgn.alignedQuery[k] != '-'){						
						qPos++;						
					}
					else{ 
						oppositeDelAlgnInfo.push_back(vector<int>{gPos, 1});//deletion information for aligning the sequences in the REVERSE ORDER (i.e. aligning the i'th sequence to the j'th sequence)

						numMut = gPos - del;
						if (numMut == 1){
							insAlgnInfo[insAlgnInfo.size() - 1][1]++; //insertion information for aligning the j'th sequence to the i'th sequence)
						}
						else{
							insAlgnInfo.push_back(vector<int>{qPos, 1}); //insertion information for aligning the j'th sequence to the i'th sequence)
						}
						del = gPos;
					}				
					currentD = qPos - gPos;
					if (currentD < minD)
						minD = currentD;
					if (currentD > maxD)
						maxD = currentD;					
				}
				clusterPairAlgnInfo[i][j].minDiag = abs(minD - clusterPairAlgnInfo[i][j].bestDiag);
				clusterPairAlgnInfo[i][j].maxDiag = abs(maxD - clusterPairAlgnInfo[i][j].bestDiag);

				clusterPairAlgnInfo[j][i].minDiag = clusterPairAlgnInfo[i][j].maxDiag;
				clusterPairAlgnInfo[j][i].maxDiag = clusterPairAlgnInfo[i][j].minDiag;

				clusterPairAlgnInfo[i][j].totalMutations = pairAlgn.numSeqIns[0] + pairAlgn.numSeqDel[0];
				clusterPairAlgnInfo[j][i].totalMutations = pairAlgn.numSeqIns[0] + pairAlgn.numSeqDel[0];
				
				clusterPairAlgnInfo[i][j].insertionEvents = insAlgnInfo;
				clusterPairAlgnInfo[j][i].insertionEvents = oppositeInsAlgnInfo;

				clusterPairAlgnInfo[i][j].deletionEvents = delAlgnInfo;
				clusterPairAlgnInfo[j][i].deletionEvents = oppositeDelAlgnInfo;

			}
			else{
				clusterPairAlgnInfo[i][j].fftMutations = 0;
				clusterPairAlgnInfo[i][j].totalMutations= 0;
				clusterPairAlgnInfo[i][j].minDiag = 0;
				clusterPairAlgnInfo[i][j].maxDiag = 0;
			}
		}		
		aligned_to_clusters[i].alignmentToOtherClusters = clusterPairAlgnInfo[i];		
	}
	

	for (int i = 0; i < numClusters; i++){
		std::free(clustersComplex[i]);
		std::free(clustersComplexQuery[i]);
		std::free(clustersFFT[i]);
		std::free(clustersFFTQuery[i]);
	}

	std::free(clustersComplexQuery);
	std::free(clustersComplex);
	std::free(clustersFFTQuery);
	std::free(clustersFFT);
	std::free(temp);
}

void QueryAlignment::ResetQueryInfo(const structvars::FASTQinfo & seqRead, int start, int end){
	/*copy sequence info to query variable*/
	numClusterHits = 0;
	num_without_gaps = 0;
	num_with_gaps = 0;
	best_cluster_so_far = NULL;
	
		
	int desiredLenToAlgn = end - start + maxClusterLength;//; querySeq[0].seq.length() + maxClusterLength, i = 0;

	currentPlan = 0;
	while (desiredLenToAlgn > fftAlgnLens[currentPlan])
		currentPlan++;

	int algnQueryLen = fftAlgnLens[currentPlan] - maxClusterLength;

	full_seq_alignment = &(*list_of_fft_plans[currentPlan][FS]);
	full_seq_alignment_RC = &(*list_of_fft_plans[currentPlan][RC]);
	
	//currentQuerySeq = querySeq[i];	
	querySeq[FS].seq = seqRead.seq.substr(start, end);
	querySeq[RC].seq = dnafunctions::ReverseComplement(querySeq[FS].seq);
	querySeq[FS].header = seqRead.seqHeader;
	querySeq[RC].header = seqRead.seqHeader;
	querySeq[FS].seqStart = start;
	querySeq[RC].seqStart = seqRead.seq.length() - end; //reverse complement
	queryLen = querySeq[FS].seq.length();
	
	/*Pad the front of the query with NOT_NT characters. All sequences willl be the same length after padding*/
	dnafunctions::PadSequence(querySeq[0].seq, algnQueryLen, FRONT);//dddddddddd
	dnafunctions::PadSequence(querySeq[1].seq, algnQueryLen, FRONT);//dddddddddddd
	
	queryStart = algnQueryLen - queryLen; //seqRead.seq.length();
	queryEnd = algnQueryLen;

	//querySeq[0].seq.insert(0, queryStart, NOT_NT);
	//querySeq[1].seq.insert(0, queryStart, NOT_NT);

	/*Create the FFT sequence variables*/
	ffthelper::convertSeq(querySeq[FS].complSeq, querySeq[FS].seq); //convert string to complex numbersdddddddddddd
	ffthelper::convertSeq(querySeq[RC].complSeq, querySeq[RC].seq); //convert string to complex numbersdddddddddd
	ffthelper::PadSeqZerosFront(querySeq[FS].complSeq, querySeq[FS].complSeqPad, algnQueryLen, fftAlgnLens[currentPlan]);
	ffthelper::PadSeqZerosFront(querySeq[RC].complSeq, querySeq[RC].complSeqPad, algnQueryLen, fftAlgnLens[currentPlan]);

	(*full_seq_alignment).ForwardFFT(querySeq[FS].complSeqPad, querySeq[FS].fourSeq); //FFT ddddddd
	(*full_seq_alignment_RC).ForwardFFT(querySeq[RC].complSeqPad, querySeq[RC].fourSeq); //FFTddddddddd

	for (int i = 0; i < maxQueryLength; i++) //reset stored positions from peptide FFTs
		seqPosStore[i] = 0;
	maxCorrectedScore = 0;
	/*for (int k = 0; k < numClusters; k++){
		for (int l = 0; l < 5; l++)
		aligned_to_clusters[k].alignToSegment[l] = true;
		}*/
}

void QueryAlignment::FindGermlineStart(){
	int alignLen = maxGermlineLength+maxClusterLength;
	FFTAlign gapplessAlignment(alignLen, algnSettings.fftParams);

	complex_t *germlinesComplex, *germlinesFFT, **clustersComplex, **clustersFFT, *temp;
	clustersComplex = (complex_t**)malloc(numClusters*sizeof(complex_t*));
	clustersFFT = (complex_t**)malloc(numClusters*sizeof(complex_t*));
	temp = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
	germlinesComplex = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
	germlinesFFT = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
	int germIndex;

	for (int j = 0; j < numClusters; j++){
		clustersComplex[j] = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
		clustersFFT[j] = (complex_t*)malloc(gapplessAlignment.alignmentLength*sizeof(complex_t));
		ffthelper::convertSeq(temp, aligned_to_clusters[j].consensus); //convert string to complex numbers	
		ffthelper::PadSeqZerosEnd(temp, clustersComplex[j], aligned_to_clusters[j].clusterSeqLength, gapplessAlignment.alignmentLength);
		gapplessAlignment.ForwardFFT(clustersComplex[j], clustersFFT[j]);		

		for (int i = 0; i < aligned_to_clusters[j].numGermlineMembers; i++){
			germIndex = aligned_to_clusters[j].germlineIndices[i];
			germlineDatabase[germIndex].seqStart = 0;
			germlineHits[germIndex].clusterIndex = j;
			
			ffthelper::convertSeq(temp, germlineDatabase[germIndex].germlineSeqData.seq); //convert string to complex numbers	
			ffthelper::PadSeqZerosEnd(temp, germlinesComplex, germlineDatabase[germIndex].germlineSeqData.seqLength, gapplessAlignment.alignmentLength);
			gapplessAlignment.ForwardFFT(germlinesComplex,germlinesFFT);
			for (int k = 0; k < gapplessAlignment.alignmentLength; k++){
				germlinesFFT[k][0] = germlinesFFT[k][0] / gapplessAlignment.alignmentLength;
				germlinesFFT[k][1] = -1 * germlinesFFT[k][1] / gapplessAlignment.alignmentLength;
			}
			gapplessAlignment.AlignSequences(clustersFFT[j],germlinesFFT); //find the best alignment diagonal
			germlineDatabase[germIndex].seqStart = gapplessAlignment.bestDiag;
			aligned_to_clusters[j].germlineStart[i] = germlineDatabase[germIndex].seqStart;
			//HERE WE CAN ACTUALLY MODEL WHERE THERE ARE MUTATIONS BETWEEEN THE GERMLINE AND CLUSTER AND STORE THEM IN ANOTHER VARIABLE
		}
	}

	for (int i = 0; i < numClusters; i++){
		std::free(clustersComplex[i]);		
		std::free(clustersFFT[i]);		
	}

	std::free(clustersComplex);	
	std::free(clustersFFT);
	std::free(germlinesComplex);
	std::free(germlinesFFT);
	std::free(temp);		
}
/*PRIVATE FUNCTIONS*/

void QueryAlignment::InitializeTransformVariables(bool allowSubstrings){		
	maxQueryLength = queryInfo.maxSeqLen;

	int minAlgnLen, avgAlgnLen = queryInfo.avgSeqLen + maxClusterLength, maxAlgnLen = maxQueryLength + maxClusterLength, meanFFTLenSize;
	if (allowSubstrings)
		minAlgnLen = minFFTAlgnScore + maxClusterLength; //OK so the fft alignment lengths will range from "maxClusterLen (assume that now query sequence below the "minFFTAlgnScore" will be allowed to pass through/have a strong alignment)" to "maxQueryLen+maxClusterLen", but on average, the majority of the lengths should be "avgSeqLen+maxClusterLen"
	else
		minAlgnLen = queryInfo.minSeqLen + maxClusterLength;
	
	printf("\tDetermining the ideal alignment size for this germline set...\n");
	fftAlgnLens.clear();
	if (queryInfo.lengthDistribution.size() == 0) //there was no length distribution calculated for this query file...then we will just make oru own bins of possible sizes without knowing the length distributino of the query file submitted
	{
		//We want to create multiple FFT models so that we can optimize the speed for alignments of different legnths. for example, we dont want a 64 base nt sequence to have the same alignment speed as a 490 base nt sequence.
		//First we should find the ideal alignment length for the "mean algn size". Most sequences will use this FFT plan. 
		
		fftAlgnLens.push_back(ffthelper::FindIdealAlgnSize(avgAlgnLen));
		meanFFTLenSize = fftAlgnLens[0];

		int i = minAlgnLen, algnLen, foundLen;
		while (true){
			algnLen = ffthelper::FindIdealAlgnSize(i);
			foundLen = false;
			for (int j = 0; j < fftAlgnLens.size(); j++){
				if (fftAlgnLens[j] == algnLen)
					foundLen = true;
			}
			if (!foundLen)
				fftAlgnLens.push_back(algnLen);
			if (i == maxAlgnLen)
				break;
			i = i*1.5;
			if (i>maxAlgnLen)
				i = maxAlgnLen;
		}
	}
	else{ //we actually determined the length distribution of the input file when reading in the sequences
		bool largestSeqPlanCreated = false;
		double perSeqs = 1,currentPercent;
		int index = 0, binLen;
		while (perSeqs > 0.1&&index < queryInfo.lengthDistribution.size()){ //lets make individual FFT plans for binned lengths which make up more than 90% of the total sequences
			binLen = queryInfo.lengthDistribution[index][2] / queryInfo.lengthDistribution[index][3]; //average length in this bin
			fftAlgnLens.push_back(ffthelper::FindIdealAlgnSize(binLen + maxClusterLength));
			currentPercent = (1.0*queryInfo.lengthDistribution[index][3]) / queryInfo.numSeqs;
			perSeqs -= currentPercent;
			index++;
		}
		
		for (int i = 0; i < fftAlgnLens.size(); i++){//lets make sure we at least accounted for hte largest sequence in the file
			if (fftAlgnLens[i] >= maxAlgnLen)
				largestSeqPlanCreated = true;
		}
		if (!largestSeqPlanCreated) //if we do not have an fft plan for that length, then lets create one
			fftAlgnLens.push_back(ffthelper::FindIdealAlgnSize(maxAlgnLen));
		meanFFTLenSize = fftAlgnLens[0];
		sort(fftAlgnLens.begin(), fftAlgnLens.end());
		
		int currLen = minAlgnLen, algnLen;		
		while (currLen<queryInfo.minSeqLen + maxClusterLength){ //ok if we allow for substrings, we have to add in extra FFT plans which account for sequence substrings < min seq length		
			algnLen = ffthelper::FindIdealAlgnSize(currLen);
			if (algnLen < fftAlgnLens[0])
				fftAlgnLens.push_back(algnLen);
			else
				break;			
			currLen = currLen*1.5;
		}
	}

	sort(fftAlgnLens.begin(), fftAlgnLens.end());
	
	/*if (!findIdealAlgnSize)
	alignmentLength = minAlignmentLength;
	else{
	alignmentLength = ffthelper::FindIdealAlgnSize(minAlignmentLength); //the given alignment length is probably not ideal speed for FFT, so we will look for a nearby length
	//alignmentLength = minAlignmentLength;
	//printf("DONT FORGET TO CHANGE ALIGNMENT LENGTH SETTING IN CREATE FFT PLAN");
	}*/
	vector<int> seqModel, germModel;
	numFFTPlans = fftAlgnLens.size();
	maxQueryLength = fftAlgnLens[numFFTPlans - 1] - maxClusterLength; //update max query length, this si because the max fft plan length might be larger than the query length
	list_of_fft_plans = new FFTAlign**[numFFTPlans];
	int s = 0, j =0;
	for (int i = 0; i < numFFTPlans; i++){		
		seqModel.resize(fftAlgnLens[i]);
		germModel.resize(fftAlgnLens[i]);
		for (int j = 0; j < maxClusterLength; j++){
			germModel[j] = j + 1;
		}
		for (int j = maxClusterLength; j < fftAlgnLens[i]; j++){
			seqModel[j] = j-maxClusterLength+1;			
		}
		
		list_of_fft_plans[i] = new FFTAlign*[2];
		list_of_fft_plans[i][FS] = new FFTAlign(fftAlgnLens[i], seqModel,germModel, algnSettings.fftParams);
		list_of_fft_plans[i][RC] = new FFTAlign(fftAlgnLens[i], seqModel, germModel, algnSettings.fftParams);//copy the plan made in the full seq alignment plan to the reverse complement plan
	}

	//full_seq_alignment = new FFTAlign(maxClusterLength + maxQueryLength, algnSettings.fftParams);
	//full_seq_alignment_RC = new FFTAlign(*full_seq_alignment); 

	int peptideFFTSize = ffthelper::FindIdealAlgnSize(algnSettings.peptideParams.peptide_len + 2 * (algnSettings.fftParams.maxGap + algnSettings.peptideParams.length_flank_gap));
	peptide_seq_alignment = new FFTAlign(peptideFFTSize, algnSettings.fftParams);
	printf("\tThe ideal alignment size for this germline set was found to be: %i \n\tA total of %i FFTW plans have been made for germline analysis\n\n", meanFFTLenSize,numFFTPlans);

	InitializeQueryVars();
	CopyQueryPrivateVarsToGermlineClusterVars();
	CreateFFTSequenceDatabase();
	//ModelOverlap();
	finalized_sw_alignments = new SWAlignment(maxQueryLength, maxClusterLength, algnSettings.swParams);
	results_no_gaps = new GermlineCluster *[numClusters];
	results_with_gaps = new GermlineCluster *[numClusters];
}

void QueryAlignment::CopyQueryPrivateVarsToGermlineClusterVars(){
	for (int i = 0; i < numClusters; i++){
		aligned_to_clusters[i].paddedNTSequence = aligned_to_clusters[i].consensus;
		//aligned_to_clusters[i].clusterEnd = aligned_to_clusters[i].consensus.length();
		dnafunctions::PadSequence(aligned_to_clusters[i].paddedNTSequence, maxClusterLength, END);
		aligned_to_clusters[i].AddSWAlignmentInfo(maxQueryLength, maxClusterLength);

		aligned_to_clusters[i].pointer_to_fft_seq_peptides = &seqPeptideList;
		//aligned_to_clusters[i].pointer_to_full_seq_alignment = &(*full_seq_alignment);
		aligned_to_clusters[i].pointer_to_peptide_alignment = &(*peptide_seq_alignment);
		aligned_to_clusters[i].pointer_to_seq_peptide_list = &seqPosStore;
		aligned_to_clusters[i].seqDel = new int*[5 * algnSettings.swParams.maxGap];
		aligned_to_clusters[i].seqIns = new int*[5 * algnSettings.swParams.maxGap];
		for (int k = 0; k < 5 * algnSettings.swParams.maxGap; k++){
			aligned_to_clusters[i].seqDel[k] = new int[2];
			aligned_to_clusters[i].seqIns[k] = new int[2];
		}
	}
}

void QueryAlignment::InitializeQueryVars(){
	
	querySeq[FS].complSeq = (complex_t*)calloc(fftAlgnLens[numFFTPlans-1]-maxClusterLength,sizeof(complex_t));
	querySeq[FS].complSeqPad = (complex_t*)calloc(fftAlgnLens[numFFTPlans-1] ,sizeof(complex_t));
	querySeq[FS].fourSeq = (complex_t*)calloc(fftAlgnLens[numFFTPlans - 1] , sizeof(complex_t));

	querySeq[RC].complSeq = (complex_t*)calloc(fftAlgnLens[numFFTPlans - 1] - maxClusterLength, sizeof(complex_t));
	querySeq[RC].complSeqPad = (complex_t*)calloc(fftAlgnLens[numFFTPlans - 1]*2 , sizeof(complex_t));
	querySeq[RC].fourSeq = (complex_t*)calloc(fftAlgnLens[numFFTPlans - 1]*2 , sizeof(complex_t));
	InitializeSeqPeptideVar();
}

void QueryAlignment::InitializeSeqPeptideVar(){
	if (seqPeptideList != NULL)
		DeleteSeqPeptideVars();
	seqPeptideList = (complex_t**)malloc(maxQueryLength*sizeof(complex_t*));
	seqPosStore = (int*)calloc(maxQueryLength, sizeof(int));
	for (int i = 0; i < maxQueryLength; i++)
		seqPeptideList[i] = (complex_t*)malloc((*peptide_seq_alignment).alignmentLength*sizeof(complex_t));
}

void QueryAlignment::DeleteSeqPeptideVars(){
	if (seqPeptideList != NULL){
		for (int i = 0; i < maxQueryLength; i++)
			std::free(seqPeptideList[i]);
		free(seqPeptideList);
	}
	seqPeptideList = NULL;

	if (seqPosStore != NULL)
		free(seqPosStore);
	seqPosStore = NULL;
}

/*
seq_peptide_data.max_gap_length = alignment.numMaxGaps; //max gaps allowed for alignment
seq_peptide_data.peptide_len = alignment.peptide_size; //size of peptides we use to chekc for gaps in different regions of sequeqnce
seq_peptide_data.length_flank_gap = alignment.extra_gap; //extra gap around gaps we allow for buffering
ffthelper::CreateFFTPeptideDatabase(seq_peptide_data, tempSeqPointer, germline_cluster.size(), alignment.germlineLength, alignment.sequenceLength); //Now make a database of FFT signals using the peptide information
alignment.swalign_data = new structvars::SWAlignment(alignment.sequenceLength, alignment.germlineLength, alignment.numMaxGaps, alignment.swParams.matchScore, alignment.swParams.mismatchScore, alignment.swParams.swGapOpen, alignment.swParams.swExtendGap);
*/
void QueryAlignment::DeleteGermlineVars(){
	if (germlineHits != NULL)
		delete[] germlineHits;
	if (finalized_sw_alignments != NULL)
		delete finalized_sw_alignments;
	if (germlineSequences != NULL)
		delete[] germlineSequences;
}

void QueryAlignment::DeleteClusterVars(){
	if (aligned_to_clusters != NULL)
		delete[] aligned_to_clusters;
}

//this function will open an input file and parse the file to update the germline cluster variable/struct
//column 1 = id/name of cluster
//column 2= all members or germline sequences that are present in this cluster
//column 3= consensus sequence
//column 4= a complex structure to define the "PairwiseClusterAlignment" structure defined in structvars.h.  This variable defines how this cluster sequence aligns to other clusters and where insertions and deletions might occur
//for this column:
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
void QueryAlignment::ReadInClusterDatabase(const string & cluster_file_name){
	std::vector<GermlineCluster> tempCluster;
	ifstream pFile;
	
	structvars::Fileinfo clusterdata;
	clusterdata.numSeqs = 0;
	GermlineCluster temp;
	std::istringstream ss;
	std::istringstream ss_cluster_members;
	string token; //use this for istringsstream, see below
	int id, current_cluster = 0;
	string cluster_members, consensus, single_member;
	vector<int> cluster_member_ids; //the actual vector of integers connecting cluster member to position in germline database file
	int minLen = 1000, maxLen = 0,totalLen =0;

	/*First quickly determine the number of clusters*/
	pFile.open(cluster_file_name.c_str());
	numClusters = 0;
	string line;
	while (!pFile.eof()){
		getline(pFile, line);
		if (line == "")
			continue;
		readfiles::RemoveTrailingSpaces(line);
		ss.clear();
		ss.str(line);
		try{
			getline(ss, token, '\t'); //first value is id
			getline(ss, token, '\t'); //first value is id
			getline(ss, token, '\t'); //first value is id
			numClusters++;
		}
		catch (exception & e){
			continue;//ignore weird sequences
		}
	}
	pFile.close();
	assert(numClusters > 0);

	aligned_to_clusters = new GermlineCluster[numClusters];

	//Each cluster will be read line by line
	pFile.open(cluster_file_name.c_str());
	while (!pFile.eof()){
		getline(pFile, line);
		readfiles::RemoveTrailingSpaces(line);
		ss.clear();
		ss.str(line);
		cluster_member_ids.clear();		
		try{
			//break down each line by tabs:
			getline(ss, token, '\t'); //first value is id
			id = stoi(token); //convert to integer
			getline(ss, token, '\t'); //next value is cluster members
			cluster_members = token;
			getline(ss, token, '\t'); //final value is consensus
			consensus = token;
			//now we need to breakdown the cluster_members by ','
			ss_cluster_members.clear();
			ss_cluster_members.str(cluster_members);
			while (getline(ss_cluster_members, single_member, ',')){
				//single_member.erase(remove_if(single_member.begin(), single_member.end(), isspace), single_member.end());
				cluster_member_ids.push_back(stoi(single_member));
			}

			if (consensus.length() < minLen)
				minLen = consensus.length();
			if (consensus.length()>maxLen)
				maxLen = consensus.length();
			totalLen += consensus.length();
			//ok, now lets parse the more complex variable structure for the 'pairwise cluster alignment' field
			int pairwiseclusterCount = 0;

			/*while (getline(ss, token, '\t')){
				tempAlParam.push_back(readfiles::ParsePairwiseAlignmentString(token));
				pairwiseclusterCount += 1;
			}*/

			aligned_to_clusters[current_cluster].consensus = consensus;
			aligned_to_clusters[current_cluster].clusterSeqLength = consensus.length();
			aligned_to_clusters[current_cluster].clusterIndex = id;
			aligned_to_clusters[current_cluster].AddGermlineMemberInfo(cluster_member_ids, germlineDatabase);
			//aligned_to_clusters[current_cluster].AddPairwiseClusterInfo(tempAlParam);
			aligned_to_clusters[current_cluster].algnParams = algnSettings;

			current_cluster++;
		}
		catch (exception & e){
			continue;//ignore weird sequences
		}
	}
	assert(current_cluster == numClusters);

	clusterSeqInfo.numSeqs = numClusters;
	clusterSeqInfo.maxSeqLen = maxLen;
	clusterSeqInfo.minSeqLen = minLen;
	clusterSeqInfo.avgSeqLen = (totalLen / numClusters) + 1;
	maxClusterLength = clusterSeqInfo.maxSeqLen; //"max germline length" just means the maximum length of the database we use to run FFT, so in fact its max cluster length
}

//If a user does not want to use clusters for the sequence data, then we can use this function to simply copy all of the data from the germline database to the cluster database
void QueryAlignment::CopyGermlineDBToClusterDB(){
	clusterSeqInfo = germlineSeqInfo;
	maxClusterLength = clusterSeqInfo.maxSeqLen;
	numClusters = numGermlines;
	aligned_to_clusters = new GermlineCluster[numGermlines];
	std::vector<int> germlineMembers(1);
	for (int i = 0; i < numGermlines; i++){
		aligned_to_clusters[i].clusterIndex = i;
		aligned_to_clusters[i].consensus = germlineDatabase[i].germlineSeqData.seq;
		aligned_to_clusters[i].clusterSeqLength = germlineDatabase[i].germlineSeqData.seq.length();
		germlineMembers[0] = i;
		aligned_to_clusters[i].AddGermlineMemberInfo(germlineMembers, germlineDatabase);
		aligned_to_clusters[i].algnParams = algnSettings;
		aligned_to_clusters[i].locusName = germlineDatabase[i].locus;
	}
	
	//printf("UH OH YOU HAVENT ADDED IN TEH FUNCTION FOR PAIRWISE ALIGNMENTS BETWEEN CLUSTESR");
	//EXIT_FAILURE;
}

//read in a file containing the germline database. add the germline sequences to the germline db vector
//FORMAT: 
	   //FIRST COLUMN => HEADER NAME OF GENE (GENE NAME, ALLELE NAME, ETC)
	   //SECOND COLUMN => ACTUAL FULL LENGTH SEQUENCE
	   //THIRD COLUMN => LOCUS OF THE GERMLINE SEQUENCE
	   //ALL OTHER COLUMNS CAN CONTAIN ADDITIONAL ANNOTATION DATA FOR GERMLINE GENE (I.E. FR1, CDR1, CHAIN, SPECIES, ETC)
void QueryAlignment::ReadGermlineFile(const string & germline_file_name, map<string, int> & unique_locus_names, map<string, int> & unique_chain_names){

	map<string, int> frameworkMap;
	for (int i = 0; i < annotationFields.size(); i++){
		frameworkMap[annotationFields[i]] = i;
	}

	ifstream pFile;
	pFile.open(germline_file_name.c_str());
	string line;
	int startpos = 0;
	int numseqs = 0;
	std::istringstream ss;
	string token; //use this for istringsstream, see below

	int minLen = 1000, maxLen = 0;

	for (int j = 0; j<germlineDatabase.size(); j++){
		if (germlineDatabase[j].seqLength < minLen)
			minLen = germlineDatabase[j].seqLength;
		if (germlineDatabase[j].seqLength>maxLen)
			maxLen = germlineDatabase[j].seqLength;			
	}

	//First lets check to see if we have a headerrow:
	string headerrow;
	getline(pFile, headerrow);
	readfiles::RemoveTrailingSpaces(headerrow);

	//break down header row into tabs
	ss.clear();
	ss.str(headerrow);
	string field;
	bool foundvheader = false;
	int seq = 0, genename = 1; //these are the regions that we will check to see if a headerline was provided
	int numVFeaturesFound = 0;
	int column = 0;
	map<int, string> annotated_vals;
	int shift;

	//we will always assume there is a header row for tab delimited data of the germline database, so this will skip the first row
	while (getline(ss, field, '\t')){//read through header row and check if we find any matches
		std::transform(field.begin(), field.end(), field.begin(), ::toupper); //convert string to uppercase

		/*
		if (field == "SEQUENCE"){
			seq = column;
			numVFeaturesFound++;
		}
		else if (field == "GENENAME")
		{
			genename = column;
			numVFeaturesFound++;
		}
		else if (field == "SHIFT"){  //we dont use this anymore
			shift = column; // SHIFT FIELD MEANS THE GERMLINE STARTS IN A DIFFERENT LOCATION THAN THE EXPECTED POSITION "1"
		}*/		
		annotated_vals[column] = field;

		column++;
	}
	/*if (numVFeaturesFound == 2)
		foundvheader = true; //all of the features we nneded are provided
	else{
		foundvheader = false; //if we did not find any features for what we are looking for, then we assume there is no header row and we just start reading in sequences
	}*/

	map<string, structvars::Abregion>::iterator my_iter;
	

	//Each cluster will be read line by line
	while (!pFile.eof()){
		getline(pFile, line);
		readfiles::RemoveTrailingSpaces(line);
		if (line != ""){
			germlineDatabase.push_back(structvars::Germline());
			numseqs = germlineDatabase.size() - 1;
			ss.clear();
			ss.str(line);
			//if (foundvheader){
			
			startpos = 0;

			getline(ss, field, '\t');//read genename
			std::transform(field.begin(), field.end(), field.begin(), ::toupper); //convert string to uppercase			
			germlineDatabase[numseqs].genename = field;//set the gene name of the germline (first column is always genename)
			getline(ss, field, '\t');//read sequence
			std::transform(field.begin(), field.end(), field.begin(), ::toupper); //convert string to uppercase			
			germlineDatabase[numseqs].germlineSeqData.seq = field;//second column of row is always sequence
			germlineDatabase[numseqs].germlineSeqData.seqLength = field.length();
			
			germlineDatabase[numseqs].locus = "";//default the locus of a germline to be empty
			germlineDatabase[numseqs].chain = "";//default the chain of a germline to be empty
			column = 2;

			while (getline(ss, field, '\t')){
				std::transform(field.begin(), field.end(), field.begin(), ::toupper); //convert string to uppercase
				string fieldname = annotated_vals[column];

				if (field != ""){
					/*if (fieldname == "SEQUENCE"){
						germlineDatabase[numseqs].germlineSeqData.seq = field;
						germlineDatabase[numseqs].germlineSeqData.seqLength = field.length();
					}
					else if (fieldname == "SHIFT") //not used anymore
						germlineDatabase[numseqs].seqStart = std::stoi(field);
					else if (fieldname == "GENENAME")
						germlineDatabase[numseqs].genename = field;
					else{*/
					if (fieldname == "FR1" || fieldname == "FR2" || fieldname == "FR3" || fieldname == "CDR1" || fieldname == "CDR2" || fieldname == "CDR3"){
						germlineDatabase[numseqs].annotation[fieldname].sequence = field;						
					}
					else if (fieldname == "LOCUS"){
						if (field != ""){
							germlineDatabase[numseqs].locus = field; //update the locus for this germline
							unique_locus_names[field]++;
						}
					}
					else if (fieldname == "CHAIN"){
						if (field != ""){
							germlineDatabase[numseqs].chain = field; //update the locus for this germline 						
							unique_chain_names[field]++;
						}
					}
					else{
						germlineDatabase[numseqs].additional_info[fieldname] = field;
					}							
					//}
				}

				column++;
			}			

			//}
			//else{
				//we assume that the first column is the sequence and the second column is the genename.
			//	getline(ss, germlineDatabase[numseqs].germlineSeqData.seq, '\t');
			//	getline(ss, germlineDatabase[numseqs].genename, '\t');
			//	std::transform(germlineDatabase[numseqs].germlineSeqData.seq.begin(), germlineDatabase[numseqs].germlineSeqData.seq.end(), germlineDatabase[numseqs].germlineSeqData.seq.begin(), ::toupper); //convert string to uppercase				
			//}
			germlineDatabase[numseqs].seqLength = germlineDatabase[numseqs].germlineSeqData.seq.length();
			germlineDatabase[numseqs].annotationIndex.resize(germlineDatabase[numseqs].seqLength); //update the index position of the germline. for 
			germlineDatabase[numseqs].additional_info["Frame"] = '1';
			int min_fr = 0;
			int max_fr = annotationFields.size()-1;
			int minp = germlineDatabase[numseqs].seqLength;
			int maxp = 0;
			std::vector<std::string> fieldnamesannotation {"FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3" };
			std::string specific_field_annotation;
			map<string, structvars::Abregion>::iterator itgerm;
			//for (auto my_iter = germlineDatabase[numseqs].annotation.begin(); my_iter != germlineDatabase[numseqs].annotation.end(); ++my_iter) {
			int find_sub_seq = 0;
			for (int gs = 0; gs < fieldnamesannotation.size(); gs++){
				specific_field_annotation = fieldnamesannotation[gs];
				itgerm = germlineDatabase[numseqs].annotation.find(specific_field_annotation);
				if (itgerm == germlineDatabase[numseqs].annotation.end())
					continue;
				//OK now that we found all of the features, we have to find there nucleotide position along the full-length germline sequence				
				readfiles::FindAbRegion(germlineDatabase[numseqs].annotation[specific_field_annotation], germlineDatabase[numseqs].germlineSeqData.seq, find_sub_seq); //SEARCH FOR FR1 IN SEQUENCE, UPDATE STARTPOS, OR POSITION AT WHICH IT WAS FOUND				

				for (int p = germlineDatabase[numseqs].annotation[specific_field_annotation].startpos; p <= germlineDatabase[numseqs].annotation[specific_field_annotation].endpos; p++){
					germlineDatabase[numseqs].annotationIndex[p] = frameworkMap[specific_field_annotation];
					if (p <= minp){
						minp = p;
						min_fr = frameworkMap[specific_field_annotation];
					}
					if (p >= maxp){
						maxp = p;
						max_fr = frameworkMap[specific_field_annotation];
					}
				}
			}

			int p = 0;
			while (p < germlineDatabase[numseqs].seqLength){
				if (germlineDatabase[numseqs].annotationIndex[p] > 0)
					break;
				germlineDatabase[numseqs].annotationIndex[p] = min_fr;
				p++;
			}

			p = germlineDatabase[numseqs].annotationIndex.size() - 1;
			while (p >= 0){
				if (germlineDatabase[numseqs].annotationIndex[p] <annotationFields.size()-1)
					break;
				germlineDatabase[numseqs].annotationIndex[p] = max_fr;
				p--;
			}
			
			germlineDatabase[numseqs].id = numseqs;
			germlineDatabase[numseqs].seqLength = germlineDatabase[numseqs].germlineSeqData.seq.length();
			if (germlineDatabase[numseqs].seqLength < minLen)
				minLen = germlineDatabase[numseqs].seqLength;
			if (germlineDatabase[numseqs].seqLength>maxLen)
				maxLen = germlineDatabase[numseqs].seqLength;			
		}
	}
	
	germlineSeqInfo.maxSeqLen = maxLen;	
	germlineSeqInfo.minSeqLen = minLen;
	germlineSeqInfo.numSeqs = germlineDatabase.size();
	numGermlines = germlineDatabase.size();
	minFFTAlgnScore = algnSettings.fftParams.scoreCutoff; //*germlineSeqInfo.minSeqLen;
	//printf("%i\n", minFFTAlgnScore);
	if (minFFTAlgnScore > round(germlineSeqInfo.maxSeqLen * 0.9)){
		printf("Warning minimum alingment score set too high based on provided germlines. Changing score to %i\n", int(round(germlineSeqInfo.maxSeqLen * 0.9)));
		minFFTAlgnScore = round(germlineSeqInfo.maxSeqLen * 0.9);
	}	
	uniqueChains.resize(unique_chain_names.size());
	uniqueLocus.resize(unique_locus_names.size());
}



//pass in: 1) a plan , the two arrays (in adn out) asociated with the input and output vectors for the plan (what the plan executes)
//pass in: 2) a list or vector of sequences
//pass in: 3) a vector of ABREAD variables which store the sequence, fourier transform, and complex variable
void QueryAlignment::CreateFFTSequenceDatabase(){
	for (int i = 0; i < numClusters; ++i){
		aligned_to_clusters[i].InitializeComplex(fftAlgnLens);
		for (int j = 0; j < numFFTPlans; j++){
			ffthelper::convertSeq(aligned_to_clusters[i].complexSeqPad, aligned_to_clusters[i].paddedNTSequence); //convert each sequence to integer
			ffthelper::PadSeqZerosEnd(aligned_to_clusters[i].complexSeqPad, aligned_to_clusters[i].complexSeqPad, maxClusterLength, (*list_of_fft_plans[j][FS]).alignmentLength); //padd integer to match the alignment len size
			(*list_of_fft_plans[j][FS]).ForwardFFT(aligned_to_clusters[i].complexSeqPad, aligned_to_clusters[i].fourierSeq[j]);
			for (int k = 0; k < (*list_of_fft_plans[j][FS]).alignmentLength; k++){
				aligned_to_clusters[i].fourierSeq[j][k][0] = aligned_to_clusters[i].fourierSeq[j][k][0] / (*list_of_fft_plans[j][FS]).alignmentLength; //(for later use with the inverse transform, fftw does not notmralize inverse fft by N
				aligned_to_clusters[i].fourierSeq[j][k][1] = -1 * aligned_to_clusters[i].fourierSeq[j][k][1] / (*list_of_fft_plans[j][FS]).alignmentLength; //(for later use with the inverse transform, fftw does not notmralize inverse fft by N...also take the conjugate (thus the negative)
			}
		}
		aligned_to_clusters[i].MakePeptideFFT(maxClusterLength);
	}
}

int QueryAlignment::FindBestAlignmentDiagonal(int guessDir, const string & locus){
	numUniqueChainHits = 0;
	numUniqueLocusHits = 0;
	numClusterHits = 0;
	num_with_gaps = 0;
	num_without_gaps = 0;
	cluster_cutoff = 0;

	int dirList[2] = { 0, 0 };
	double inverseRatio = (1 / algnSettings.scoreFoldRatio), scoreW;
	
	double maxCS = 0, maxCL = 0;

	bool findMultiplePeaksF, findMultiplePeaksR, foundDir = false;
	FFTAlign *winningStrand;
	int minAlgnLen;
	//	double otherHitData[7];
	//int minAlignmentPosition, maxAlignmentPosition, numSubElements;
	clock_t s1 = clock();
	for (int j = 0; j < numClusters; j++){
		minAlgnLen = min(queryLen, aligned_to_clusters[j].clusterSeqLength);
		findMultiplePeaksF = true;
		findMultiplePeaksR = true;

		if (locus != "" && aligned_to_clusters[j].locusName!="" && aligned_to_clusters[j].locusName!=locus){//if parameter defining the locus is supplied, then only consider clusters whose locus matches the input parameter
			aligned_to_clusters[j].FFTgappedScore = -100;
			aligned_to_clusters[j].FFTgaplessScore = -100;
			aligned_to_clusters[j].SWgapPenalty = 0;

			continue;
		}

		if (!foundDir  && dirList[0]>algnSettings.numObsAboveFoldRatio && dirList[1] == 0){
			guessDir = 1;//if you see that the first three times resulted in the forward direction, then guess direction = forward always
			foundDir = true;
		}
		else if (!foundDir && dirList[1] > algnSettings.numObsAboveFoldRatio && dirList[0] == 0){
			guessDir = -1; //if you see that the first three times resulted in the RC direction, then guess direction = reverseend
			foundDir = true;
		}

		//####RUN FFT ALIGNMENT#######/
		if (guessDir == 0 || guessDir == 1)
			(*full_seq_alignment).AlignSequences(querySeq[0].fourSeq, aligned_to_clusters[j].fourierSeq[currentPlan]);
		else{
			(*full_seq_alignment).bestDiagScore = 0.01;
			(*full_seq_alignment).gappedScore = 0.0;
			findMultiplePeaksF = false;
		}

		if (guessDir == 0 || guessDir == -1)
			(*full_seq_alignment_RC).AlignSequences(querySeq[1].fourSeq, aligned_to_clusters[j].fourierSeq[currentPlan]);
		else{
			(*full_seq_alignment_RC).bestDiagScore = 0.01;
			(*full_seq_alignment_RC).gappedScore = 0.0;
			findMultiplePeaksR = false;
		}
		//#############################//

		scoreW = (*full_seq_alignment).bestDiagScore / (*full_seq_alignment_RC).bestDiagScore;

		if (scoreW < 1){// fold score comparing best forward alignment to RC alignment//if less than 1, then the RC had a better score
			if (scoreW <= inverseRatio){
				dirList[1]++;
				(*full_seq_alignment).bestDiagScore = 0.01;
				(*full_seq_alignment).gappedScore = 0.0;
				findMultiplePeaksF = false;
			}
		}
		else{
			if (scoreW >= algnSettings.scoreFoldRatio){
				dirList[0]++;
				(*full_seq_alignment_RC).bestDiagScore = 0.01;
				(*full_seq_alignment_RC).gappedScore = 0.0;
				findMultiplePeaksR = false;
			}
		}

		if (findMultiplePeaksF)
			(*full_seq_alignment).IdentifyAdditionalAlignments(minAlgnLen);

		if (findMultiplePeaksR)
			(*full_seq_alignment_RC).IdentifyAdditionalAlignments(minAlgnLen);

		if ((*full_seq_alignment).gappedScore >= (*full_seq_alignment_RC).gappedScore){
			winningStrand = full_seq_alignment;
			aligned_to_clusters[j].dir = FS;// '+';
			//aligned_to_clusters[j].queryInfoPointer = &querySeq[aligned_to_clusters[j].dir];// &querySeq[0];
		}
		else{
			winningStrand = full_seq_alignment_RC;
			aligned_to_clusters[j].dir = RC; //'-';
			//aligned_to_clusters[j].queryInfoPointer = &querySeq[1];
		}
		aligned_to_clusters[j].queryInfoPointer = &querySeq[aligned_to_clusters[j].dir];
		aligned_to_clusters[j].pointer_to_full_seq_alignment = winningStrand;
		aligned_to_clusters[j].querySeqStart = queryStart;
		aligned_to_clusters[j].querySeqEnd = queryEnd;
		aligned_to_clusters[j].minDiag = (*winningStrand).minDiag;
		aligned_to_clusters[j].maxDiag = (*winningStrand).maxDiag;
		aligned_to_clusters[j].numFFTGaps = (*winningStrand).numGaps;
		aligned_to_clusters[j].containsGaps = aligned_to_clusters[j].numFFTGaps > 0;
		aligned_to_clusters[j].FFTgaplessScore = (*winningStrand).bestDiagScore;
		aligned_to_clusters[j].FFTgappedScore = (*winningStrand).gappedScore;
		aligned_to_clusters[j].bestDiag = (*winningStrand).bestDiag;
		aligned_to_clusters[j].numSteps = 0; //delete me
		aligned_to_clusters[j].SWgapPenalty = 0;
		if (aligned_to_clusters[j].FFTgappedScore > maxCS){
			maxCS = aligned_to_clusters[j].FFTgappedScore;
			maxCL = (*winningStrand).alignmentModel[(*winningStrand).bestDiag][5];
		}
		winningStrand = NULL;
	}

	if (maxCS < minFFTAlgnScore){
	//if (maxCS <  algnSettings.fftParams.scoreCutoff*maxCL){
		numClusterHits = 0;
		return 0;
	}
	clock_t s2 = clock();
	t += s1 - s2;
	cluster_cutoff = max(4.0, round(algnSettings.fftParams.cluster_threshold*maxCS));
	double bestAlignmentNoGap = 0;
	int swap_no_gap_index = 0;

	//organize the winning clusters or those above threshold into clusters that contain gaps and clusters taht do not contain gaps
	for (int i = 0; i < numClusters; i++){
		if (aligned_to_clusters[i].FFTgappedScore >= cluster_cutoff){
			numClusterHits++;
			aligned_to_clusters[i].UpdateAlignmentToSequence((*full_seq_alignment).alignmentModel[aligned_to_clusters[i].bestDiag]);
			if (aligned_to_clusters[i].containsGaps){
				num_with_gaps++;
				results_with_gaps[num_with_gaps - 1] = &aligned_to_clusters[i];
			}
			else{
				num_without_gaps++;
				results_no_gaps[num_without_gaps - 1] = &aligned_to_clusters[i];
				if (aligned_to_clusters[i].FFTgaplessScore > bestAlignmentNoGap){
					bestAlignmentNoGap = aligned_to_clusters[i].FFTgaplessScore;
					swap_no_gap_index = num_without_gaps - 1;
				}
			}
		}
	}

	if (num_without_gaps > 0){
		best_cluster_so_far = results_no_gaps[swap_no_gap_index];
		results_no_gaps[swap_no_gap_index] = results_no_gaps[0];
		results_no_gaps[0] = best_cluster_so_far;
	}
	else{
		best_cluster_so_far = NULL;
	}

	return numClusterHits;
}

void QueryAlignment::OptimizeVGeneClusterAlignmentWithGaps(int seqAnalyzed){
	clock_t a, b;
	a = clock();
	if (num_without_gaps > 0){

		(*best_cluster_so_far).FindAllGapsInQueryToClusterAlgn(); //make sure there are no gaps in this sequence
		
		if ((*best_cluster_so_far).finalScore > maxCorrectedScore)
			maxCorrectedScore = (*best_cluster_so_far).finalScore; //update the best alignment score we currently can find

		for (int j = 1; j < num_without_gaps; j++){
			if (results_no_gaps[j]->dir == best_cluster_so_far->dir && results_no_gaps[j]->bestDiag >= best_cluster_so_far->bestDiag - algnSettings.fftParams.maxGap && results_no_gaps[j]->bestDiag <= best_cluster_so_far->bestDiag + algnSettings.fftParams.maxGap)
				(*results_no_gaps[j]).FindAllGapsUsingAnotherRefCluster(*best_cluster_so_far);
			else
				(*results_no_gaps[j]).FindAllGapsInQueryToClusterAlgn();
			if ((*results_no_gaps[j]).finalScore > maxCorrectedScore)
				maxCorrectedScore = (*results_no_gaps[j]).finalScore; //update the best alignment score we currently can find
		}

		for (int j = 0; j < num_with_gaps; j++){
			if (results_with_gaps[j]->dir == best_cluster_so_far->dir && results_with_gaps[j]->bestDiag >= best_cluster_so_far->bestDiag - algnSettings.fftParams.maxGap && results_with_gaps[j]->bestDiag <= best_cluster_so_far->bestDiag + algnSettings.fftParams.maxGap)
				(*results_with_gaps[j]).FindAllGapsUsingAnotherRefCluster(*best_cluster_so_far);
			else
				(*results_with_gaps[j]).FindAllGapsInQueryToClusterAlgn();
			if ((*results_with_gaps[j]).finalScore > maxCorrectedScore)
				maxCorrectedScore = (*results_with_gaps[j]).finalScore; //update the best alignment score we currently can find
		}
	}
	else{
		for (int j = 0; j < num_with_gaps; j++){
			(*results_with_gaps[j]).FindAllGapsInQueryToClusterAlgn();
			if ((*results_with_gaps[j]).finalScore > maxCorrectedScore)
				maxCorrectedScore = (*results_with_gaps[j]).finalScore; //update the best alignment score we currently can find
		}
	}

	//DebugAlignments(seqAnalyzed);

	b = clock();
	tO += (b - a);
	//now select the clusters that were found with 70% of the max cluster//
	cluster_cutoff = max(0.7*minFFTAlgnScore, 0.7*maxCorrectedScore);
	cluster_cutoff = max(15.0, cluster_cutoff);
}

void QueryAlignment::DebugAlignments(int seqNum){
	if (seqNum == 1){
		debugfile << "SeqNum" << "\t" << "FFTGAPSINCLUSTER" << "\t" << "loopnumber" << "\t" << "clusternumber" << "\t" << "Total Insertions from algorithm" << "\t" << "Total Deletions from algorithm" << "\t" << "Gap Penalty from algorithm" << "\t" << "max score from algorithm" << "\t" << "Total Insertions from full alignment" << "\t" << "Total Deletions from full alignment" <<"\t" << "gap penalty from fulll alignment" << "\t" << "max score from full alignment" << "\t" << "gapless query from algorithm" << "\t" << "gapless query from full alignment" << "\n";
	}

	for (int i = 0; i < num_without_gaps; i++){
		(*finalized_sw_alignments).OverlapAlignComplete(results_no_gaps[i]->queryInfoPointer->seq,results_no_gaps[i]->paddedNTSequence);
		debugfile << seqNum << "\t" << "NOGAPS" << "\t" << i << "\t" << results_no_gaps[i]->clusterIndex << "\t" << results_no_gaps[i]->numSeqIns[0] << "\t" << results_no_gaps[i]->numSeqDel[0] << "\t" << results_no_gaps[i]->SWgapPenalty << "\t" << results_no_gaps[i]->finalScore << "\t" << finalized_sw_alignments->numSeqIns[0] << "\t" << finalized_sw_alignments->numSeqDel[0] << "\t" << finalized_sw_alignments->gapPenalties << "\t" << finalized_sw_alignments->maxScore << "\t" << results_no_gaps[i]->querySeqToBeModified << "\t" << finalized_sw_alignments->gaplessQuery << "\n";
	}
	for (int i = 0; i < num_with_gaps; i++){
		(*finalized_sw_alignments).OverlapAlignComplete(results_with_gaps[i]->queryInfoPointer->seq, results_with_gaps[i]->paddedNTSequence);
		debugfile << seqNum << "\t" << "WITHGAPS"<<"\t"<< i << "\t" << results_with_gaps[i]->clusterIndex << "\t" << results_with_gaps[i]->numSeqIns[0]<<"\t"<<results_with_gaps[i]->numSeqDel[0] << "\t" <<results_with_gaps[i]->SWgapPenalty<<"\t"<< results_with_gaps[i]->finalScore << "\t" <<finalized_sw_alignments->numSeqIns[0]<<"\t"<<finalized_sw_alignments->numSeqDel[0]<<"\t"<< finalized_sw_alignments->gapPenalties << "\t" << finalized_sw_alignments->maxScore << "\t" << results_with_gaps[i]->querySeqToBeModified << "\t" << finalized_sw_alignments->gaplessQuery << "\n";
	}
}

void QueryAlignment::OptimizeJGeneClusterAlignmentWithGaps(){
	clock_t a, b;
	a = clock();
	for (int j = 0; j < num_without_gaps; j++){
		//(*results_no_gaps[j]).FindAllGapsInQueryToClusterAlgn();
		(*results_no_gaps[j]).SearchGapsInFrontandEnd();
		if ((*results_no_gaps[j]).finalScore > maxCorrectedScore)
			maxCorrectedScore = (*results_no_gaps[j]).finalScore; //update the best alignment score we currently can find
	}
	for (int j = 0; j < num_with_gaps; j++){
		(*results_with_gaps[j]).SearchGapsInFrontandEnd(); //FindAllGapsInQueryToClusterAlgn();
		if ((*results_with_gaps[j]).finalScore > maxCorrectedScore)
			maxCorrectedScore = (*results_with_gaps[j]).finalScore; //update the best alignment score we currently can find
	}
			
	b = clock();
	tO += (b - a);
	//now select the clusters that were found with 90% of the max cluster//
	cluster_cutoff = max(0.7*minFFTAlgnScore,0.7*maxCorrectedScore);
	cluster_cutoff = max(4.0, cluster_cutoff);
}

int QueryAlignment::AlignClusterResultsToGermline(){
	string seq_to_align, gene_seq;
	double gapPenalty;
	
	int foundHits = 0;
	//	double minVal;
	double scoreBest;

	
	int algnDiagData[4];
	//double *topHits = new double[algnSettings.maxHits]();
	double minScore = 0;
	int shift;
	clock_t a, b;
	a = clock();
	int temp;
	int locusCount, chainCount;
	bool newLocus, newChain;
	for (int i = 0; i < num_without_gaps; i++){
		UpdateUniqueLocus(results_no_gaps[i]->locusName);				
		if (results_no_gaps[i]->finalScore >= cluster_cutoff){
			seq_to_align = results_no_gaps[i]->querySeqToBeModified;
			gapPenalty = results_no_gaps[i]->SWgapPenalty;
			for (int j = 0; j < results_no_gaps[i]->numGermlineMembers; j++){
				//if (results_no_gaps[i]->germlineIndices[j] == 198)
					//cout << "stop";
				shift = results_no_gaps[i]->germlineStart[j];
				germlineHits[results_no_gaps[i]->germlineIndices[j]].swAlignDiag = results_no_gaps[i]->initialSWDiag + shift;
				algnDiagData[0] = results_no_gaps[i]->alignmentToSequence[1];
				algnDiagData[1] = results_no_gaps[i]->alignmentToSequence[3];
  				if (shift != 0){
					shift -= algnDiagData[1];
					if (shift >= 0){
						algnDiagData[1] = 0;
						algnDiagData[0] += shift;
					}
					else{
						algnDiagData[1] = -1 * shift;
					}
				}

				shift = algnDiagData[0] - (algnDiagData[1] - 0);
					
				algnDiagData[1] = shift < queryStart ? queryStart - shift : 0; //make this always above next line (algnDiagdata'0] = shift..)!!!
				algnDiagData[0] = shift < queryStart ? queryStart : shift;																			
				algnDiagData[2] = results_no_gaps[i]->germlineSequences[j].length()-algnDiagData[1];//results_no_gaps[i]->alignmentToSequence[5];


					
				(*finalized_sw_alignments).GaplessAlignDiag(seq_to_align, results_no_gaps[i]->germlineSequences[j], algnDiagData);
				scoreBest = (*finalized_sw_alignments).maxScore;
				scoreBest -= gapPenalty;

				if (scoreBest < cluster_cutoff)
					//This is mostlikely due to a the fact that this germline is in the cluster gruop BUT with respect to query it is truncated.
					//so the part of the cluster that aligns to query IS NOT found in this germline
					continue;

				if ((*finalized_sw_alignments).germStart < 10 && (*finalized_sw_alignments).germStart>0)
				{
					shift = (*finalized_sw_alignments).seqStart - (*finalized_sw_alignments).germStart;
					(*finalized_sw_alignments).germStart = shift < queryStart ? queryStart - shift: 0;
					(*finalized_sw_alignments).seqStart = shift < queryStart ? queryStart : shift;
				}
						
				shift = (results_no_gaps[i]->germlineSequences[j].length()-1) - (*finalized_sw_alignments).germEnd;
				if (shift > 0 && shift < 10){						
				
					(*finalized_sw_alignments).germEnd = shift + (*finalized_sw_alignments).seqEnd >= querySeq[0].seq.length() ? (*finalized_sw_alignments).germEnd + ((querySeq[0].seq.length()-1) - (*finalized_sw_alignments).seqEnd) : results_no_gaps[i]->germlineSequences[j].length()-1;
					(*finalized_sw_alignments).seqEnd = shift + (*finalized_sw_alignments).seqEnd >= querySeq[0].seq.length() ? querySeq[0].seq.length() - 1 : shift + (*finalized_sw_alignments).seqEnd;
				}
					
				aligned_to_germline_scores[foundHits][0] = scoreBest;
				aligned_to_germline_scores[foundHits][1] = results_no_gaps[i]->germlineIndices[j];					
				germlineHits[results_no_gaps[i]->germlineIndices[j]].score = scoreBest;
				germlineHits[results_no_gaps[i]->germlineIndices[j]].minDiag = results_no_gaps[i]->minDiag;
				germlineHits[results_no_gaps[i]->germlineIndices[j]].maxDiag = results_no_gaps[i]->maxDiag;
				germlineHits[results_no_gaps[i]->germlineIndices[j]].dir = results_no_gaps[i]->dir;
				germlineHits[results_no_gaps[i]->germlineIndices[j]].startGerm = (*finalized_sw_alignments).germStart;
				germlineHits[results_no_gaps[i]->germlineIndices[j]].endGerm = (*finalized_sw_alignments).germEnd;
					
				germlineHits[results_no_gaps[i]->germlineIndices[j]].startQuery = (*finalized_sw_alignments).seqStart;
				germlineHits[results_no_gaps[i]->germlineIndices[j]].endQuery = (*finalized_sw_alignments).seqEnd;

				germlineHits[results_no_gaps[i]->germlineIndices[j]].algnLen = germlineHits[results_no_gaps[i]->germlineIndices[j]].endQuery - germlineHits[results_no_gaps[i]->germlineIndices[j]].startQuery + 1;
					
				foundHits++;
			}
		}
	}

	for (int i = 0; i < num_with_gaps; i++){
		UpdateUniqueLocus(results_with_gaps[i]->locusName);
		if (results_with_gaps[i]->finalScore >= cluster_cutoff){
			seq_to_align = results_with_gaps[i]->querySeqToBeModified;
			gapPenalty = results_with_gaps[i]->SWgapPenalty;
			
				for (int j = 0; j < results_with_gaps[i]->numGermlineMembers; j++){
					shift = results_with_gaps[i]->germlineStart[j];
					germlineHits[results_with_gaps[i]->germlineIndices[j]].swAlignDiag = results_with_gaps[i]->initialSWDiag + shift;
					algnDiagData[0] = results_with_gaps[i]->alignmentToSequence[1];
					algnDiagData[1] = results_with_gaps[i]->alignmentToSequence[3];
					if (shift != 0){
						shift -= algnDiagData[1];
						if (shift >= 0){
							algnDiagData[1] = 0;
							algnDiagData[0] += shift;
						}
						else{
							algnDiagData[1] = -1 * shift;
						}
					}

					shift = algnDiagData[0] - (algnDiagData[1] - 0);

					algnDiagData[1] = shift < queryStart ? queryStart - shift : 0; //make this always above next line (algnDiagdata'0] = shift..)!!!

					algnDiagData[0] = shift < queryStart ? queryStart : shift;


					algnDiagData[2] = results_with_gaps[i]->germlineSequences[j].length();//results_no_gaps[i]->alignmentToSequence[5];

					
					(*finalized_sw_alignments).GaplessAlignDiag(seq_to_align, results_with_gaps[i]->germlineSequences[j], algnDiagData);
					scoreBest = (*finalized_sw_alignments).maxScore;
					scoreBest -= gapPenalty;

					if ((*finalized_sw_alignments).germStart < 10 && (*finalized_sw_alignments).germStart>0)
					{
						shift = (*finalized_sw_alignments).seqStart - (*finalized_sw_alignments).germStart;
	
						(*finalized_sw_alignments).germStart = shift < queryStart ? queryStart-shift : 0;
						(*finalized_sw_alignments).seqStart = shift < queryStart ? queryStart : shift;
					}

					shift = (results_with_gaps[i]->germlineSequences[j].length()-1) - (*finalized_sw_alignments).germEnd;
					if (shift > 0 && shift < 10){
						
						(*finalized_sw_alignments).germEnd = shift + (*finalized_sw_alignments).seqEnd >= querySeq[0].seq.length() ? (*finalized_sw_alignments).germEnd + ((querySeq[0].seq.length() - 1) - (*finalized_sw_alignments).seqEnd) : results_with_gaps[i]->germlineSequences[j].length()-1;
						(*finalized_sw_alignments).seqEnd = shift + (*finalized_sw_alignments).seqEnd >= querySeq[0].seq.length() ? querySeq[0].seq.length() - 1 : shift + (*finalized_sw_alignments).seqEnd;
					}

					aligned_to_germline_scores[foundHits][0] = scoreBest;
					aligned_to_germline_scores[foundHits][1] = results_with_gaps[i]->germlineIndices[j];					
					germlineHits[results_with_gaps[i]->germlineIndices[j]].score = scoreBest;
					germlineHits[results_with_gaps[i]->germlineIndices[j]].dir = results_with_gaps[i]->dir;
					germlineHits[results_with_gaps[i]->germlineIndices[j]].startGerm = (*finalized_sw_alignments).germStart;
					germlineHits[results_with_gaps[i]->germlineIndices[j]].endGerm = (*finalized_sw_alignments).germEnd;

					germlineHits[results_with_gaps[i]->germlineIndices[j]].minDiag = results_with_gaps[i]->minDiag;
					germlineHits[results_with_gaps[i]->germlineIndices[j]].maxDiag = results_with_gaps[i]->maxDiag;

					germlineHits[results_with_gaps[i]->germlineIndices[j]].startQuery = (*finalized_sw_alignments).seqStart;
					germlineHits[results_with_gaps[i]->germlineIndices[j]].endQuery = (*finalized_sw_alignments).seqEnd;
					germlineHits[results_with_gaps[i]->germlineIndices[j]].algnLen = germlineHits[results_with_gaps[i]->germlineIndices[j]].endQuery - germlineHits[results_with_gaps[i]->germlineIndices[j]].startQuery + 1;

					foundHits++;
				}
		}
	}
	
	if (foundHits == 0)
		return 0;
	
	sort(aligned_to_germline_scores.begin(), aligned_to_germline_scores.begin() + foundHits , [](const vector<double> & a, const vector<double> & b){ return (a[0] > b[0]); });
	int num_unique_scores = 1;
	int k = 1;

		
	while (num_unique_scores<=algnSettings.maxHits && k < foundHits){		
		if (aligned_to_germline_scores[k][0] < aligned_to_germline_scores[k - 1][0])
			num_unique_scores++;		
		k++;
	}

	numGermlineHits = num_unique_scores > algnSettings.maxHits ? k-1 : k;

	AlignGermlineToQuery();

	b = clock();
	tO += (b - a);
	return numGermlineHits;
}

void QueryAlignment::AlignGermlineToQuery(){
	int Diagcoordinates[7];
	int index;
	string query,germline;
	
	
	for (int i = 0; i < numGermlineHits; i++){ //for the top top germlines found to have a good alignment to the query, perform a proper alignment to the query
		index = aligned_to_germline_scores[i][1];
		
		query = aligned_to_clusters[germlineHits[index].clusterIndex].queryInfoPointer->seq;				
		germline = germlineSequences[index];

		Diagcoordinates[0] = germlineHits[index].startGerm;
		Diagcoordinates[1] = germlineHits[index].endGerm;
		Diagcoordinates[2] = germlineHits[index].swAlignDiag;//we need to use the INITIAL DIAGONAL, NOT THE MODIFIED ONE. this is because we are alignment to the query, not the gap-modified query //germlineHits[index].startQuery - germlineHits[index].startGerm;
		Diagcoordinates[3] = germlineHits[index].minDiag;
		Diagcoordinates[4] = germlineHits[index].maxDiag;		
		Diagcoordinates[5] = queryStart;
		Diagcoordinates[6] = germline.length();

		(*finalized_sw_alignments).OverlapAlignDiag(query, germline,Diagcoordinates);
		(*finalized_sw_alignments).TracebackWithAnnotation();
		
		germlineHits[index].numDel = (*finalized_sw_alignments).numSeqDel[0];
		germlineHits[index].numIns = (*finalized_sw_alignments).numSeqIns[0];

		germlineHits[index].numMatch = (*finalized_sw_alignments).totalMatch;
		germlineHits[index].numMismatch = (*finalized_sw_alignments).totalMismatch;
		germlineHits[index].startQuery = (*finalized_sw_alignments).seqStart;
		germlineHits[index].endQuery = (*finalized_sw_alignments).seqEnd;
		germlineHits[index].startGerm = (*finalized_sw_alignments).germStart;
		germlineHits[index].endGerm = (*finalized_sw_alignments).germEnd;
		germlineHits[index].algnLen = germlineHits[index].endQuery - germlineHits[index].startQuery + 1;
		germlineHits[index].germlineAlSeq = (*finalized_sw_alignments).alignedGerm;
		germlineHits[index].queryAlSeq = (*finalized_sw_alignments).alignedQuery;
		germlineHits[index].score = (*finalized_sw_alignments).maxScore;
		
		AnnotateQuery(index, query);

		aligned_to_germline_scores[i][0] = germlineHits[index].score;		
	}
	sort(aligned_to_germline_scores.begin(), aligned_to_germline_scores.begin() + numGermlineHits, [](const vector<double> & a, const vector<double> & b){ return (a[0] > b[0]); });
	int num_unique_scores = 1;
	int k = 1;

	while (num_unique_scores <= algnSettings.maxHits && k < numGermlineHits){
		if (aligned_to_germline_scores[k][0] < aligned_to_germline_scores[k - 1][0])
			num_unique_scores++;
		k++;
	}

	numGermlineHits = num_unique_scores > algnSettings.maxHits ? k - 1 : k;

}

void QueryAlignment::AnnotateQuery(int index, const string & query){	
	int start, end;
	bool nextRegion = true;
	
	germlineHits[index].algnRegion[0] = germlineDatabase[index].annotationIndex[(*finalized_sw_alignments).germStart] > 0 && germlineDatabase[index].annotationIndex[(*finalized_sw_alignments).germStart]  <7 ? germlineDatabase[index].annotationIndex[(*finalized_sw_alignments).germStart] : 0;
	germlineHits[index].algnRegion[1] = germlineDatabase[index].annotationIndex[(*finalized_sw_alignments).germEnd] > 0 && germlineDatabase[index].annotationIndex[(*finalized_sw_alignments).germEnd] < 7 ? germlineDatabase[index].annotationIndex[(*finalized_sw_alignments).germEnd] : 7;
	
	for (int i = 0; i < 8; i++){
		germlineHits[index].included[annotationFields[i]] = false;
	}	
	
	for (int i = germlineHits[index].algnRegion[0]; i <= germlineHits[index].algnRegion[1]; i++){				
		if (germlineDatabase[index].annotation.find(annotationFields[i]) != germlineDatabase[index].annotation.end()){
			germlineHits[index].included[annotationFields[i]] = true;
			start = germlineDatabase[index].annotation[annotationFields[i]].startpos >= (*finalized_sw_alignments).germStart ? germlineDatabase[index].annotation[annotationFields[i]].startpos : (*finalized_sw_alignments).germStart;
			end = germlineDatabase[index].annotation[annotationFields[i]].endpos <= (*finalized_sw_alignments).germEnd ? germlineDatabase[index].annotation[annotationFields[i]].endpos : (*finalized_sw_alignments).germEnd;
			
			germlineHits[index].annotation[annotationFields[i]].startpos = (*finalized_sw_alignments).tracebackData[start][SEQPOS];// : (*finalized_sw_alignments).seqStart; //(*finalized_sw_alignments).tracebackData[(*finalized_sw_alignments).germStart][SEQPOS];
			germlineHits[index].annotation[annotationFields[i]].endpos = (*finalized_sw_alignments).tracebackData[end][SEQPOS];//(*finalized_sw_alignments).tracebackData[germlineDatabase[index].annotation[annotationFields[i]].endpos][SEQPOS] : (*finalized_sw_alignments).seqEnd;// tracebackData[(*finalized_sw_alignments).germStart][SEQPOS];// (*finalized_sw_alignments).tracebackData[germlineDatabase[index].annotation[annotationFields[i]].endpos][SEQPOS];
			germlineHits[index].annotation[annotationFields[i]].startpos_algnseq = (*finalized_sw_alignments).tracebackData[start][ALGNPOS];// : (*finalized_sw_alignments).seqStart; //(*finalized_sw_alignments).tracebackData[(*finalized_sw_alignments).germStart][SEQPOS];
			germlineHits[index].annotation[annotationFields[i]].endpos_algnseq = (*finalized_sw_alignments).tracebackData[end][ALGNPOS];//(*finalized_sw_alignments).tracebackData[germlineDatabase[index].annotation[annotationFields[i]].endpos][SEQPOS] : (*finalized_sw_alignments).seqEnd;// tracebackData[(*finalized_sw_alignments).germStart][SEQPOS];// (*finalized_sw_alignments).tracebackData[germlineDatabase[index].annotation[annotationFields[i]].endpos][SEQPOS];
			germlineHits[index].annotation[annotationFields[i]].sequence = query.substr(germlineHits[index].annotation[annotationFields[i]].startpos, germlineHits[index].annotation[annotationFields[i]].endpos - germlineHits[index].annotation[annotationFields[i]].startpos+1);
			germlineHits[index].match[annotationFields[i]] = (*finalized_sw_alignments).tracebackData[end + 1][MATCHES] - (*finalized_sw_alignments).tracebackData[start][MATCHES];
			germlineHits[index].mismatch[annotationFields[i]] = (*finalized_sw_alignments).tracebackData[end + 1][MISMATCHES] - (*finalized_sw_alignments).tracebackData[start][MISMATCHES];
			germlineHits[index].indel[annotationFields[i]] = (*finalized_sw_alignments).tracebackData[end + 1][INDEL] - (*finalized_sw_alignments).tracebackData[start][INDEL];
		}
	}
}

void QueryAlignment::SummarizeVGeneResults(structvars::VGeneResults & results){
	int index = int(aligned_to_germline_scores[0][1]);
	results.direction = bestSeqDir[germlineHits[index].dir];
	results.QS = germlineHits[index].startQuery - queryStart + querySeq[germlineHits[index].dir].seqStart;
	results.QE = germlineHits[index].endQuery - queryStart + querySeq[germlineHits[index].dir].seqStart;
	results.GS = germlineHits[index].startGerm;
	results.GE = germlineHits[index].endGerm;
	results.queryAlSeq = germlineHits[index].queryAlSeq;
	results.germlineAlSeq = germlineHits[index].germlineAlSeq;
	results.germlineCodingFrame = std::stoi(germlineDatabase[index].additional_info["Frame"]);
	results.totalMatch = germlineHits[index].numMatch;
	results.totalMismatch = germlineHits[index].numMismatch;
	results.totalIndel = germlineHits[index].numDel + germlineHits[index].numIns;
	
	string geneNames = "", geneScores = "";

	results.startingAnnotation = annotationFields[germlineDatabase[index].annotationIndex[results.GS]];
	results.endingAnnotation = annotationFields[germlineDatabase[index].annotationIndex[results.GE]];
	results.codonStart = dnafunctions::AdjustToCorrectFrame(results.QS, results.GS, results.germlineCodingFrame);

	bool first_of=false;
	int start,end;
	
	for (int i = 1; i < 7; i++){
		if (germlineHits[index].included[annotationFields[i]]){
			
			results.sequence[annotationFields[i]] = germlineHits[index].annotation[annotationFields[i]].sequence;
			if (annotationFields[i] == "CDR3"){
				results.cdr3start = germlineHits[index].annotation[annotationFields[i]].startpos - queryStart + querySeq[germlineHits[index].dir].seqStart;
				results.startPos[annotationFields[i]] = results.cdr3start;// std::to_string(results.cdr3start);
				results.startPosText[annotationFields[i]] = std::to_string(results.cdr3start);
				if (results.startingAnnotation == "CDR3")
					start = results.cdr3start;
			}
			else{
				start = germlineHits[index].annotation[annotationFields[i]].startpos - queryStart + querySeq[germlineHits[index].dir].seqStart;
				results.startPos[annotationFields[i]] = start;
				results.startPosText[annotationFields[i]] = std::to_string(start);
			}

			results.startPos_algnseq[annotationFields[i]] = std::to_string(germlineHits[index].annotation[annotationFields[i]].startpos_algnseq);
			results.endPos_algnseq[annotationFields[i]] = std::to_string(germlineHits[index].annotation[annotationFields[i]].endpos_algnseq);
			end = germlineHits[index].annotation[annotationFields[i]].endpos - queryStart + querySeq[germlineHits[index].dir].seqStart;
			results.endPos[annotationFields[i]] = std::to_string(end);
			
			results.numMatch[annotationFields[i]] = std::to_string(germlineHits[index].match[annotationFields[i]]);
			results.numMisMatch[annotationFields[i]] = std::to_string(germlineHits[index].mismatch[annotationFields[i]]);
			results.numIndel[annotationFields[i]] = std::to_string(germlineHits[index].indel[annotationFields[i]]);
			results.readingFrame[annotationFields[i]] = start % 3 + 1;
			results.readingFrameText[annotationFields[i]] = std::to_string(results.readingFrame[annotationFields[i]]);
		}
		else{			
			results.sequence[annotationFields[i]] = "";
			results.endPos[annotationFields[i]] = "";
			results.startPosText[annotationFields[i]] = "";
			results.startPos[annotationFields[i]] = -1;
			results.numMatch[annotationFields[i]] = "";
			results.numMisMatch[annotationFields[i]] = "";
			results.numIndel[annotationFields[i]] = ""; 
			results.readingFrameText[annotationFields[i]] = "";
			results.readingFrame[annotationFields[i]] = -1;
			results.startPos_algnseq[annotationFields[i]] = "";
			results.endPos_algnseq[annotationFields[i]] = "";
		}		
	}

	results.readingFrame[results.startingAnnotation] = results.codonStart % 3 + 1;
	results.readingFrameText[results.startingAnnotation] = std::to_string(results.readingFrame[results.startingAnnotation]);

	for (int i = 0; i < numGermlineHits - 1; i++){
		index = int(aligned_to_germline_scores[i][1]);
		geneNames += germlineDatabase[int(aligned_to_germline_scores[i][1])].genename + ",";
		geneScores += std::to_string(aligned_to_germline_scores[i][0]) + ",";		
	}
	geneNames += germlineDatabase[int(aligned_to_germline_scores[numGermlineHits - 1][1])].genename;
	geneScores += std::to_string(aligned_to_germline_scores[numGermlineHits - 1][0]);
	results.genes = geneNames;
	results.scores = geneScores;
	
}

void QueryAlignment::SummarizeJGeneResults(structvars::JGeneResults & results){
	int index = int(aligned_to_germline_scores[0][1]);

	results.direction = bestSeqDir[germlineHits[int(aligned_to_germline_scores[0][1])].dir];
	results.QS = germlineHits[int(aligned_to_germline_scores[0][1])].startQuery - queryStart + querySeq[germlineHits[int(aligned_to_germline_scores[0][1])].dir].seqStart;
	results.QE = germlineHits[int(aligned_to_germline_scores[0][1])].endQuery - queryStart + querySeq[germlineHits[int(aligned_to_germline_scores[0][1])].dir].seqStart;
	results.GS = germlineHits[int(aligned_to_germline_scores[0][1])].startGerm;
	results.GE = germlineHits[int(aligned_to_germline_scores[0][1])].endGerm;
	results.germlineCodingFrame = std::stoi(germlineDatabase[index].additional_info["Frame"]);
	results.codonStart = dnafunctions::AdjustToCorrectFrame(results.QS, results.GS, results.germlineCodingFrame);
	results.queryAlSeq = germlineHits[index].queryAlSeq;
	results.germlineAlSeq = germlineHits[index].germlineAlSeq;
	results.totalMatch = germlineHits[index].numMatch;
	results.totalMismatch = germlineHits[index].numMismatch;
	results.totalIndel = germlineHits[index].numDel + germlineHits[index].numIns;
	string geneNames = "", geneScores = "";
	for (int i = 0; i < numGermlineHits - 1; i++){
		geneNames += germlineDatabase[int(aligned_to_germline_scores[i][1])].genename + ",";
		geneScores += std::to_string(aligned_to_germline_scores[i][0]) + ",";
	}
	geneNames += germlineDatabase[int(aligned_to_germline_scores[numGermlineHits - 1][1])].genename;
	geneScores += std::to_string(aligned_to_germline_scores[numGermlineHits - 1][0]);
	results.genes = geneNames;
	results.scores = geneScores;
}

clock_t QueryAlignment::ReturnTime(){
	return t;
}

clock_t QueryAlignment::ReturnTimeOp(){
	return tO;
}

void QueryAlignment::PrintClusterResults(const std::string & filename){
	ofstream outfile(filename.c_str());
	for (int i = 0; i < numClusters; i++){
		outfile << aligned_to_clusters[i].clusterIndex << "\t";
		for (int j = 0; j < aligned_to_clusters[i].numGermlineMembers-1; j++){
			outfile << aligned_to_clusters[i].germlineIndices[j] << ",";
		}
		outfile << aligned_to_clusters[i].germlineIndices[aligned_to_clusters[i].numGermlineMembers-1] << "\t";
		outfile << aligned_to_clusters[i].consensus << "\n";		
	}
	outfile.close();
}

void QueryAlignment::AddUniqueLociToVectorString(std::vector<std::string> & locivector, int & currentLen){ //if you pass in a vector of strings, then it updates this vector with unique Loci
	bool newLocus = true;
	for (int i = 0; i < numUniqueLocusHits; i++){
		newLocus = true;
		for (int j = 0; j < currentLen; j++){
			if (uniqueLocus[i] == locivector[j]){
				newLocus = false;
				break;
			}
		}
		if (newLocus && uniqueLocus[i] != ""){
			locivector[currentLen] = uniqueLocus[i];// .push_back(uniqueLocus[i]);
			currentLen++;
		}
	}	
}