#include <iostream>
#include "string.h"
#include <iomanip>
#include "SWAlignment.h"
#include "StructVars.h"

using namespace std;
int main(){
	structvars::SWAlignSettings clusterFormationSettings;
	clusterFormationSettings.matchScore = 5;
	clusterFormationSettings.mismatchScore = 4;
	clusterFormationSettings.swGapOpen = 10;//in order to ensure the germlines clustered together lack gaps, we want this to be very sensitive to gap openings
	clusterFormationSettings.swExtendGap = 5;

	string s2 = "TAATTTGAAGGAATTGCAAGGACGAGGA";
	string s1 = "TACTTTTGAAAATTGCAGGATGTTAAGA";
	SWAlignment pairAlgn(305, 305, clusterFormationSettings);
	int algnCoord[5];
	algnCoord[0] = 0;
	algnCoord[1] = s1.length();
	algnCoord[2] = 0;
	algnCoord[3] = s1.length();
	algnCoord[4] = s1.length();
	algnCoord[5] = 0;
	algnCoord[6] = s1.length();
	pairAlgn.OverlapAlignDiag(s2, s1, algnCoord);
	pairAlgn.TracebackWithAnnotation();
	pairAlgn.printAlignment();
	int germpos;
	while (true){
		cout << endl << endl << "Enter germline index: " << endl;
		cin >> germpos;
		if (germpos < 0)
			break;
		pairAlgn.printTracebackAnnotation(germpos);
		cout << endl << endl;
		pairAlgn.printAlignment();
	}
	cout << "program complete, press enter" << endl;
	cin.get();
	//for (int i = 0; i<1000; i++){		
	//	pairAlgn.OverlapAlignComplete(s1, s2); //align each cluster to the "seed cluster consensus". 
	//}
	return 0;
}