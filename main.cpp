#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <numeric>

using namespace std;

//struct for storing on target data
struct iden
{
	char nt1;
	char nt2;
	int position;
	double odds_score;
};

//helper functions for splitting lines on a delimeter
vector<string> split(string s, char delimiter)
{
	size_t pos = 0;
	string token;
	vector <string> splitStr;
	while ((pos = s.find(delimiter)) != string::npos)
	{
		token = s.substr(0, pos);
		splitStr.push_back(token);
		s.erase(0, pos + 1);
	}
	splitStr.push_back(s);
	return splitStr;
}

//read in input sequences
void read(string &filename, vector<vector<string>> &refSeqs, vector<vector<string>> &querySeqs)
{
	ifstream fin;
	string line;
	int cnt = 0;

	fin.open(filename);

	if (fin.is_open())
	{
		while (getline(fin, line))
		{
			if (line.find(">") != string::npos)
			{
				getline(fin, line);
				querySeqs[cnt].push_back(line);
				getline(fin, line);
				cnt++;
			}
			else
			{
				refSeqs[cnt].push_back(line);
			}
		}
	}
	fin.close();
}

//write out results of off target analysis
void write(string &filename, vector<vector<string>> &refSeqs, vector<vector<string>> &querySeqs, vector<vector<double>> &shScores, vector<vector<double>> &offTargetScores, vector<double> &queryOffTargetScores)
{
	ofstream fout;

	fout.open(filename);

	if (fout.is_open())
	{
		for (int i = 0; i < refSeqs.size(); i++)
		{
			for (int j = 0; j < refSeqs[i].size(); j++)
			{
				fout << refSeqs[i][j] << " : " << shScores[i][j] << "," << offTargetScores[i][j] << endl;
			}

			fout << ">" << endl;
			fout << querySeqs[i][0] << " : " << queryOffTargetScores[i] << endl;
			fout << ">" << endl;
		}
	}

	fout.close();
}

//load on target data from CASPERinfo
void loadCRISPRSCAN(string &filename, vector<iden> &idens)
{
	ifstream fin;
	string line;

	fin.open(filename);

	if (fin.is_open())
	{
		while (getline(fin, line))
		{
			if (line.find("DATA:CRISPRSCAN") != string::npos)
			{
				break;
			}
		}

		while (getline(fin, line))
		{
			if (line.find("------------------------------------------------------") != string::npos)
			{
				break;
			}
			else
			{
				iden nid;

				vector<string> mytoke = split(line, '\t');
				nid.nt1 = mytoke[0][0];
				nid.nt2 = mytoke[0][1];
				nid.position = stoi(mytoke[1]);
				nid.odds_score = stod(mytoke[2]);

				idens.push_back(nid);
			}
		}
	}

}

//load hsu matrix with default values
void loadHSUMatrix(map<string, vector<double>> &hsuMatrix)
{
	hsuMatrix["GT"] = { 1, 1.6533,0.9030,1.5977,0.9235,0.8070,0.9632,1.0163,0.2658,0.7119,1.5211,0.6942,1.0434,0.5255,0.8981,0.7164,0.8399,0.5376,0.2821,0.6898 };
	hsuMatrix["AC"] = { 1, 1.5142,1.1597,1.6582,0.9924,0.0247,0.5522,1.8687,0.7737,0.9270,0.7292,0.4842,0.4824,0.7060,1.0221,0.0181,0.3496,0.1811,0.1362,0.2700 };
	hsuMatrix["GG"] = { 1, 1.3234,1.4157,1.2967,1.2060,0.9524,0.2304,1.0163,0.8100,1.1559,0.7075,1.5791,0.3490,0.0899,0.0497,0.0045,0.2267,0.2153,0.5250,0.4965 };
	hsuMatrix["TG"] = { 1, 1.5366,1.2771,1.2689,1.2197,1.0645,0.7791,1.2445,0.9885,1.5319,0.7184,1.7381,0.4166,0.1285,0.0720,0.0549,0.2261,0.3119,0.1343,0.0601 };
	hsuMatrix["TT"] = { 1, 1.7347,0.8215,1.1579,1.1816,0.7380,0.9004,0.8368,0.2997,0.6210,1.1400,0.3561,0.6192,0.1799,0.2665,0.2793,0.2613,0.1152,0.1680,0.3372 };
	hsuMatrix["CA"] = { 1, 1.0186,0.9649,2.2504,1.6222,0.2405,0.7561,1.0651,0.1102,1.4293,0.3533,0.6178,0.7269,0.1165,0.0367,0.5013,0.4147,0.1786,0.5315,0.1664 };
	hsuMatrix["CT"] = { 1, 1.6719,0.9688,1.0732,1.0869,0.6475,1.0142,0.8635,0.3059,0.4487,0.9046,0.4327,0.5576,0.1379,0.0722,0.3279,0.2420,0.0433,0.1351,0.4403 };
	hsuMatrix["GA"] = { 1, 1.1662,0.4544,2.7867,1.0461,0.6036,0.8132,0.7875,0.6882,1.3655,0.1240,0.1953,0.2497,0.0132,0.0227,0.0478,0.3682,0.3175,0.5621,0.4588 };
	hsuMatrix["AA"] = { 1, 1.1916,1.0954,2.8709,1.1310,0.5160,0.6439,1.0322,0.5356,1.2868,0.0780,0.2592,0.2675,0.0469,0.0252,0.0052,0.0218,0.1718,0.6970,0.2720 };
	hsuMatrix["AG"] = { 1, 1.4786,1.0820,1.2952,0.7450,0.9763,0.4912,0.9272,0.6022,1.0375,0.3047,0.8210,0.0303,0.0365,0.0592,0.0253,0.1553,0.1006,0.2175,0.0275 };
	hsuMatrix["TC"] = { 1, 1.0400,0.9954,1.6466,1.3410,0.0102,0.5428,2.3401,0.4367,0.2143,0.3405,0.2640,0.0935,0.0462,0.0688,0.0165,0.3659,0.0546,0.0857,0.2543 };
	hsuMatrix["CC"] = { 1, 1.0345,1.0478,1.0507,1.4075,0.0540,0.6396,2.0810,0.4585,0.1555,0.1369,0.1026,0.0417,0.0105,0.0458,0.0099,0.2114,0.0552,0.0253,0.0596 };
}

//calculate on target scores for sequence
int calculateOnTargetScore(string &sequence, vector<iden> &idens)
{
	double totalScore = 0.0;
	int returnScore = 0;

	for (int i = 0; i < idens.size(); i++)
	{
		char nucleo1 = idens.at(i).nt1;  //may have to change this to pointers b/c of multiple creations of object
		char nucleo2 = idens.at(i).nt2;
		int pos = idens.at(i).position;
		if (pos < sequence.size())
		{
			if (nucleo2 != 'x')
			{
				string dinucleo = string() + nucleo1 + nucleo2;
				if (sequence.substr(pos - 1, 2) == dinucleo)
				{
					totalScore += idens.at(i).odds_score;
				}
			}
			else
			{
				if (sequence.at(pos - 1) == nucleo1)
				{
					totalScore += idens.at(i).odds_score;
				}
			}
		}
	}

	//following line normalizes to best possible score
	totalScore = 1 - ((1.29401 - totalScore) / 1.94947);
	returnScore = (totalScore * 100) + 0.5; //0.5 for proper rounding
	return returnScore;
}

// OffTarget scoring functions 
double shScore(vector<int> &mismatches, vector<string> &hsuKeys, map<string, vector<double>> &hsuMatrix)
{
	double tot_sh = 1.0;
	for (int i = 0; i < mismatches.size(); i++)
	{
		tot_sh *= hsuMatrix[hsuKeys[i]][mismatches[i]];
	}
	return tot_sh;
}

double ssScore(vector<int> &mismatches)
{
	double tot_ss = 1.0;
	for (int i = 0; i < mismatches.size(); i++)
	{
		if (mismatches[i] < 6)
		{
			tot_ss -= 0.1;
		}
		else if (mismatches[i] < 12)
		{
			tot_ss -= 0.05;
		}
		else
		{
			tot_ss -= 0.0125;
		}
	}
	return tot_ss;
}

double stScore(vector<int> &mismatches)
{
	double tot_st = 3.5477;
	for (int i = 0; i < mismatches.size(); i++)
	{
		tot_st -= 1.0 / (mismatches[i] + 1);
	}
	return tot_st / 3.5477;
}

char reverseComp(char &c)
{
	char n = 'N';
	switch (c)
	{
	case 'A':
		n = 'T';
		break;
	case 'T':
		n = 'A';
		break;
	case 'G':
		n = 'C';
		break;
	case 'C':
		n = 'G';
		break;
	}
	return n;
}

bool getMismatches(string &refSeq, string &currentQuerySeq, vector<int> &mismatchLocations, vector<string> &hsuKeys, int &seqLength, int &maxMismatches)
{
	/* vars */
	string hsuKey;

	//character by character comparison of ref and query sequences
	for (int j = seqLength - 1; j >= 0; j--)
	{
		try
		{
			if (refSeq.at(j) != currentQuerySeq.at(j))
			{
				//store mismatch location
				mismatchLocations.push_back(j);

				//store key for HSU matrix
				hsuKey = string() + currentQuerySeq.at(j) + reverseComp(refSeq.at(j));
				hsuKeys.push_back(hsuKey);
			}
		}
		catch (int e)
		{
			cout << "mismatch error" << endl;
		}

		/* if there are too many mismatches, break */
		if (mismatchLocations.size() > maxMismatches)
		{
			return false;
		}
	}
	return true;
}

double findSimilars(vector<string> &refSeqs, vector<int> &refOnScores, string &currentQuerySeq, int &currentQueryScore, vector<double> &targetScores, vector<double> &shScores, int seqLength, map<string, vector<double>> &hsuMatrix, int &maxMismatches)
{
	/* vars */
	double rRatio = 0.0;
	double value = 0.0;
	string refSeq;

	/* loop through each organims unique seq from CSPR file and compare against the given query sequence */
	for (unsigned long i = 0; i < refSeqs.size(); i++)
	{
		vector<int> mismatches;
		vector<string> hsuKeys;
		rRatio = refOnScores[i] / currentQueryScore;
		refSeq = refSeqs[i];

		//character by character comparison of ref and query sequences
		if (getMismatches(refSeq, currentQuerySeq, mismatches, hsuKeys, seqLength, maxMismatches))
		{
			//update runnning score if mismatch count wasn't too large
			value = (((sqrt(shScore(mismatches, hsuKeys, hsuMatrix)) + stScore(mismatches)) * pow(ssScore(mismatches), 6) * pow(rRatio, 2)) / 4);
			shScores.push_back(shScore(mismatches, hsuKeys, hsuMatrix));
			targetScores.push_back(value);
		}
	}


	//avg the score
	double averageScore = accumulate(targetScores.begin(), targetScores.end(), 0.0);
	if (targetScores.size() != 0)
	{
		averageScore /= (targetScores.size());
	}
	return averageScore;

}

//calculate off target score
void offTargetScore(vector<vector<string>> &refSeqs, vector<vector<int>> &refOnScores, vector<vector<string>> &querySeqs, vector<vector<int>> &queryOnScores, vector<vector<double>> &offTargetScores, vector<vector<double>> &shScores, map<string, vector<double>> &hsuMatrix, int &maxMismatches, vector<double> &queryOffTargetScores)
{
	//loop through each set of ref seqs
	for (int i = 0; i < refSeqs.size(); i++)
	{
		double OTScore = findSimilars(refSeqs[i], refOnScores[i], querySeqs[i][0], queryOnScores[i][0], offTargetScores[i], shScores[i], 20, hsuMatrix, maxMismatches);
		queryOffTargetScores.push_back(OTScore);
	}
}

int main()
{
	//vars
	string inputFilename = "C:\\Users\\Tfry\\Desktop\\CASPER_alg_seqs.txt";
	string outputFilename = "C:\\Users\\Tfry\\Desktop\\CASPER-alg-results.txt";
	string CASPERinfoFilename = "C:\\Users\\Tfry\\Desktop\\CASPERapp\\CASPERinfo";
	
	
	vector<vector<string>> refSeqs (3);
	vector<vector<string>> querySeqs (3);
	vector<vector<int>> refOnScores (3);
	vector<vector<int>> queryOnScores (3);
	vector<vector<double>> offTargetScores (3);
	vector<vector<double>> shScores (3);
	vector<iden> idens;
	map<string, vector<double>> hsuMatrix;
	vector<double> queryOffTargetScores;
	
	int maxMismatches = 5;

	//read input data
	read(inputFilename, refSeqs, querySeqs);

	//parse CASPERinfo for on target data
	loadCRISPRSCAN(CASPERinfoFilename, idens);

	//load hsu matrix
	loadHSUMatrix(hsuMatrix);

	//calculate on target scores for ref sequences
	for (int i = 0; i < refSeqs.size(); i++)
	{
		for (int j = 0; j < refSeqs[i].size(); j++)
		{
			refOnScores[i].push_back(calculateOnTargetScore(refSeqs[i][j], idens));
		}
	}

	//calculate on target scores for query sequences
	for (int i = 0; i < querySeqs.size(); i++)
	{
		for (int j = 0; j < querySeqs[i].size(); j++)
		{
			queryOnScores[i].push_back(calculateOnTargetScore(querySeqs[i][j], idens));
		}
	}

	//calculate off target scores for query sequences
	offTargetScore(refSeqs, refOnScores, querySeqs, queryOnScores, offTargetScores, shScores, hsuMatrix, maxMismatches, queryOffTargetScores);

	//write results
	write(outputFilename, refSeqs, querySeqs, shScores, offTargetScores, queryOffTargetScores);

	system("pause");
	return 0;
}