//        _
//       | | _
//     _ | || |
//    | || || |
//    | || || |


#include "clustr.h"

using namespace std;

clustr::clustr()
{
	Sigma = 0;
	CoG = 0;
	Charge = 0;
}

void clustr::SetParameters()
{
	double chargesum = 0;
	double position = 0;
	double positionSoN = 0;
	double sigmapos = 0;
	double SoNsum = 0;
	double noisesum = 0;
	
	for (auto const& i : Elements)
	{
		chargesum += i.second;
		position += i.second*i.first;
		SoNsum += (i.second/Noise.at(i.first));
		noisesum += pow(Noise.at(i.first), 2);
		positionSoN += i.first*(i.second/Noise.at(i.first));
	}
	Charge = chargesum;
	CoG = position/chargesum;
	CoGSoN = positionSoN/SoNsum;
	Significance = SoNsum;
	Significance2 = chargesum/sqrt(noisesum);
	for (auto const& i : Elements)
	{
		sigmapos +=  pow((CoG - i.first), 2)*i.second;
	}
	Sigma = sqrt(sigmapos/chargesum);
}


int clustr::MaxCharge()
{
	double charge_max = 0;
	int strip_max;
	for (auto const& i : Elements)
	{
		if (i.second > charge_max)
		{
			charge_max = i.second;
			strip_max = i.first;
		}
	}
	return strip_max;
}

int clustr::MaxSoN()
{
	double SoN_max = 0;
	int strip_max;
	//cout << "Debug MAXSON2: 1 " << endl;
	for (auto const& i : Elements)
	{
		
		double SoN = i.second/Noise.at(i.first);
		if (SoN > SoN_max)
		{
			SoN_max = SoN;
			strip_max = i.first;
		}
	}
	//cout << "Debug MAXSON2: 2 " << endl;
	return strip_max;
}

void clustr::Erase(int key)
{
	Noise.erase(key);
	Elements.erase(key);
}

clustr::~clustr() { }
