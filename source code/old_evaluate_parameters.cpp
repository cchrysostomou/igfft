//This function will evaluate the parameters inputed from the user
void EvaluateParameters(int argc, char *argv[]){
	char var;
	int currentVar = 2;
	//	int parameterVarInt;
	FILE *pFile;
	string filename; //for reading filenames/strings
	std::istringstream ss;
	string token; //use this for istringsstream, see below
	int numDatabase;
	if (argc < 4){//check input arguments
		printf("Error: You must enter the location of the sequence file, and the location of at least one database to compare to\n");
		exit(EXIT_FAILURE);
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
	variables_used.fileinputName = argv[1];//input file name
	variables_used.fileoutputName = variables_used.fileinputName + ".fft.annotation";
	while (currentVar < argc){
		var = argv[currentVar][1];
		if (strlen(argv[currentVar]) != 2 && argv[currentVar][0] != '-'){
			printf("Incorrect parameter\n");
			exit(EXIT_FAILURE);
		}
		cout << argv[currentVar] << endl;
		switch (var)
		{
		case 'v': //location of the vgermline database
			currentVar++;
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
						variables_used.fileGermlineV[numDatabase] = token; //the first index MUST be the location of the list of germline sequences.  the second index MUST be the location of the clustered germline seuqences.
						numDatabase++;
					}
				}
				catch (std::exception&e){
					printf("The V Gene Database must be defined by two file locations: the first location is the list of germline sequences, and the second file refers to the clustered germline sequences");
					exit(EXIT_FAILURE);
				}
			}
			//error out if we did not read 2 files (index starts at 0)
			if (numDatabase < 1 || variables_used.fileGermlineV[0] == "")
			{
				printf("The V Gene Database file location must be defined");
				exit(EXIT_FAILURE);
			}
			//define whether user desires vgene data
			variables_used.containsVGermline = true;
			break;
		case 'j': //location of the jgermline database
			//define whether the user desires j gene data
			variables_used.containsJGermline = true;
			currentVar++;
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
					variables_used.fileGermlineJ[numDatabase] = token; //the first index MUST be the location of the list of germline sequences.  the second index MUST be the location of the clustered germline seuqences.
					numDatabase++;
				}
				catch (std::exception&e){
					printf("The J Gene Database must be defined by two file locations: the first location is the list of germline sequences, and the second file refers to the clustered germline sequences");
					exit(EXIT_FAILURE);
				}
			}
			//error out if we did not read 2 files (index starts at 0)
			if (numDatabase < 1 || variables_used.fileGermlineJ[0] == "")
			{
				printf("The J Gene Database file location must be provided");
				exit(EXIT_FAILURE);
			}
			break;
		case 'o': //location of the file output name
			currentVar++;
			variables_used.fileoutputName = string(argv[currentVar]);
			break;
		case 'm': //match score
			currentVar++;
			try{
				variables_used.query_algn_settings.swParams.matchScore = std::stoi(argv[currentVar]);
			}
			catch (std::exception& e)
			{
				printf("Match score must be integer\n");
				exit(EXIT_FAILURE);
			}
			break;
		case 's': //mismatch score
			currentVar++;
			try{
				variables_used.query_algn_settings.swParams.mismatchScore = std::stoi(argv[currentVar]);
			}
			catch (std::exception& e)
			{
				printf("Match score must be integer\n");
				exit(EXIT_FAILURE);
			}
			break;
		case 'g': //maximum number of gaps allowed
			currentVar++;
			try{
				variables_used.max_gap = std::stoi(argv[currentVar]);
			}
			catch (std::exception& e)
			{
				printf("Gap must be integer\n");
				exit(EXIT_FAILURE);
			}
			break;
		case 'e': //extra space for aligning small peptides
			currentVar++;
			try{
				variables_used.extra_gap = std::stoi(argv[currentVar]);
			}
			catch (std::exception& e)
			{
				printf("Extra gap value must be integer\n");
				exit(EXIT_FAILURE);
			}
			break;
		case 'a': //minimum alignment length
			currentVar++;
			try{
				variables_used.minAlgnLength = std::stoi(argv[currentVar]);
			}
			catch (std::exception& e)
			{
				printf("Extra gap value must be integer\n");
				exit(EXIT_FAILURE);
			}
			break;
		case 'p': //percent identityt
			currentVar++;
			try{
				variables_used.percentIdentity = std::stod(argv[currentVar]);
			}
			catch (std::exception& e)
			{
				printf("percent identity must be a number between 0 and 1\n");
				exit(EXIT_FAILURE);
			}
			if (variables_used.percentIdentity > 1 || variables_used.percentIdentity < 0)
			{
				printf("percent identity must be a number between 0 and 1\n");
				exit(EXIT_FAILURE);
			}
			cout << variables_used.percentIdentity << endl;
			break;
		case 'c': //score cutoff
			currentVar++;
			try{
				variables_used.query_algn_settings.fftParams.scoreCutoff = std::stod(argv[currentVar]);
				if (variables_used.percentIdentity > 1 || variables_used.percentIdentity < 0)
				{
					printf("cutoff must be a number between 0 and 1\n");
					exit(EXIT_FAILURE);
				}
			}
			catch (std::exception& e)
			{
				printf("cutoff must be a number between 0 and 1\n");
				exit(EXIT_FAILURE);
			}
			break;
		case 'i': //input file format
			currentVar++;
			try{
				variables_used.query_algn_settings.inputFileFormat = string(argv[currentVar]);
				for (int k = 0; k < variables_used.query_algn_settings.inputFileFormat.length(); k++){
					variables_used.query_algn_settings.inputFileFormat[k] = ::toupper(variables_used.query_algn_settings.inputFileFormat[k]);
				}
				cout << variables_used.query_algn_settings.inputFileFormat << endl;
				if (variables_used.query_algn_settings.inputFileFormat != "FASTA" && variables_used.query_algn_settings.inputFileFormat != "TAB"){
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
			variables_used.query_algn_settings.ignoreHeader = string(argv[currentVar]) == "false" ? false : true;
			break;
		default:
			currentVar++;
		}
		currentVar++;
	}
	//check if either of the databases were included. at least one must be defined.
	//that is are they both false....
	if (!(variables_used.containsJGermline || variables_used.containsVGermline)){
		printf("Error: You must enter the location of at least one germline database\n");
		exit(EXIT_FAILURE);
	}
}