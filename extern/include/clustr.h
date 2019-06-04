#ifndef __CLUSTR_H__
#define __CLUSTR_H__
#include <iostream>
#include <string>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <iterator>
#include <ctime>

using namespace std;

class clustr
{
	public:
	// Constructor
		clustr();	
	// Deconstructor
		~clustr();
	
	
		double CoG ;
		double Charge;
		double SoN;
		double CoGSoN;
		double SigmaCoG;
		
		std::map<int,double> Elements;
		std::multimap<int,double> ElementsNoise;
		
		void SetParameters();
		void SetParameters_w_Noise();
		void SetParameters_w_Noise2();
		int MaxCharge();
		int MaxCharge_w_Noise();
		int MaxSoN();
		//void setEventList(vector<pair<double, double>>);
} ;



#endif
