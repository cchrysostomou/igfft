#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <set>
#include "DnaFunctions.h"

using namespace std;
using namespace dnafunctions;

int main(){
	
	//string aa = TranslateSeq("GNNGTGACGTTGGACGAGTCCGGGGGCGGCCTCCAGACGCCCGGAGGAGGGCTCAGCCTCGTCTGCAAGGGCTCCGGGTTCACCTTCAGCAGTTACGGCATGTACTGGGTGCGACAGGCGCCCGGCAAGGGGCTGGAATTCGTCGCGGGTATTAGCGGTGAGGGGAGTAGCACAGGATACGGGCCGGCGGTGGAGGGGCGTGCCACCATCTCGGGGGGGGACGGGGAGAGGACAGGGGGGCTGCAGCTGG",2);	
	cout << "ok";
	string aa = TranslateSeq("GGGAC",1);
	cout << aa << endl;
	cout << "done";
	cin.get();
	return 0;
	
}