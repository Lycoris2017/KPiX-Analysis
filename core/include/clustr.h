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
		double Significance; //linear addition of charge/noise for each strip
		double Significance2; //Adding noises in quadrature and later dividing chargesum by this noise sum
		double CoGSoN;
		double Sigma;
		
		std::map<int,double> Elements;
		std::map<int,double> Noise;
		
		void SetParameters();
		int MaxCharge();
		int MaxSoN();
		
		void Erase(int);
		//void setEventList(vector<pair<double, double>>);
} ;



#endif
