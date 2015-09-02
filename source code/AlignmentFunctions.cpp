

#include <iostream>

#include <vector>
#include "AlignmentFunctions.h"


using namespace std;


//model 1= sequence model
//model 2 = germline model
/*%this functino assumes => NO ALIASING BETWEEN TO SIGNALS...THAT WAY WE KNOW
%ANY OVERLAP IS ONLY DUE TO A SPECIFIC ALIGNMENT.  Any signals that overlap
%with one another would result in us not really knowing whcih overlapping sequenes
%are truly leading to the correct overlap
%i.e. compare these signal pairs [0,0,0,1,2,3];[1,2,3,0,0,0]
%WITH [1,2,3];[1,2,3];

%model1 => this defines how the nucleotide sequence that is not shifting is
%currently aligned.  it lists all the positions where a dna sequence
%exists.  For example if this sequence is put first as [ACTG.....] WHERE
%"." represents no sequence or "0 padding" then model1 would be
%[1,2,3,4,0,0,0,0,0] where 4 represents position 4 or nucleotide 4 along
%sequence

%model 2=> this defines how the nucleotide sequence that is shifting along
%the first one is aligned/designed.  For example if model 2 is [.....ACTG]
%where '." represents regions of no sequence or "0" padding, then we right
%it as [0,0,0,0,0,1,2,3,4]
*/
vector<vector<int> > algnfunct::ModelOverlap(vector<int> model1, vector<int> model2, int sizeModel)
{
	vector<int> model2shift(sizeModel);
	vector<vector<int> >algnShift(sizeModel);
	int shiftPos;	
	vector<int> algn(sizeModel);	
	int currentPos = 1;
	int start = -1;
	int maxLen = 0;
	int maxP = -1;
	int endP;
	int batch = -1;
	vector<int> temp(5);
	vector<vector<int> > posS;
	vector<int> temp2(3);

	int skip;
	int algnIndLeft, algnIndRight;
	for (int shiftI = 0; shiftI<sizeModel; ++shiftI){
		for (int i = 0; i<sizeModel; ++i){
			shiftPos = ((i - shiftI) + sizeModel) % sizeModel;
			model2shift[i] = model2[shiftPos];
		}
		for (int i = 0; i<sizeModel; ++i){
			algn[i] = model1[i] * model2shift[i];
		}
		start = -1;
		endP = -1;
		maxP = -1;
		maxLen = 0;
		skip = 0;
		for (int s = 0; s<sizeModel; ++s){
			if (algn[s] != 0){
				if (start == -1)
					start = s;
			}
			else{
				if (start != -1){
					endP = s - 1;
					batch = batch + 1;
					temp2[0] = start;
					temp2[1] = endP;
					temp2[2] = endP - start + 1;
					posS.push_back(temp2);
					if (temp2[2]>maxLen){
						maxP = batch;
						maxLen = temp2[2];
					}
					start = -1;
					endP = -1;
				}
			}
		}
		if (algn[sizeModel - 1] != 0){
			if (start == -1){
				start = sizeModel - 1;
				endP = sizeModel - 1;
			}
			else{
				endP = sizeModel - 1;
			}

			batch = batch + 1;
			temp2[0] = start;
			temp2[1] = endP;
			temp2[2] = endP - start;
			posS.push_back(temp2);
			if (temp2[2]>maxLen)
				maxP = batch;
		}
		temp[0] = shiftI;
		if (maxP == -1){
			temp[1] = 0;
			temp[2] = 0;
			temp[3] = 0;
			temp[4] = 0;
		}
		else{
			algnIndLeft = posS[maxP][0];

			algnIndRight = posS[maxP][1];
			temp[1] = model1[algnIndLeft];
			temp[2] = model1[algnIndRight];
			temp[3] = model2shift[algnIndLeft];
			temp[4] = model2shift[algnIndRight];
		}
		algnShift[shiftI] = temp;
	}
	return algnShift;
}