//        _
//       | | _
//     _ | || |
//    | || || |
//    | || || |


#include "clustr.h"

using namespace std;

clustr::clustr()
{
	SigmaCoG = 0;
	CoG = 0;
	Charge = 0;
}

void clustr::SetParameters()
{
	double chargesum = 0;
	double position = 0;
	double sigmapos = 0;
	
	for (auto const& i : Elements)
	{
		chargesum += i.second;
		position += i.second*i.first;
	}
	Charge = chargesum;
	CoG = position/chargesum;
	
	for (auto const& i : Elements)
	{
		sigmapos +=  pow((CoG - i.first), 2)*i.second;
	}
	SigmaCoG = sqrt(sigmapos/chargesum);;
}

void clustr::SetParameters_w_Noise()
{
	double chargesum = 0;
	double position = 0;
	double sigmapos = 0;
	pair<multimap<int,double>::iterator, multimap<int,double>::iterator> it_intern;
	
	for (multimap<int,double>::iterator it = ElementsNoise.begin(); it != ElementsNoise.end(); it++)
	{
		it_intern=ElementsNoise.equal_range((*it).first);
		multimap<int,double>::iterator Charge_it = it_intern.first;
		chargesum += Charge_it->second;
		position += Charge_it->second * Charge_it->first;
	}
	Charge = chargesum;
	CoG = position/chargesum;
	
	for (multimap<int,double>::iterator it = ElementsNoise.begin(); it != ElementsNoise.end(); it++)
	{
		it_intern=ElementsNoise.equal_range((*it).first);
		multimap<int,double>::iterator Charge_it = it_intern.first;

		sigmapos +=  pow((CoG - Charge_it->first), 2)*Charge_it->second;
	}
	SigmaCoG = sqrt(sigmapos/chargesum);;
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

int clustr::MaxCharge_w_Noise()
{
	double charge_max = 0;
	int strip_max;
	pair<multimap<int,double>::iterator, multimap<int,double>::iterator> it_intern;
	for (multimap<int,double>::iterator it = ElementsNoise.begin(); it != ElementsNoise.end(); it++)
	{
		it_intern=ElementsNoise.equal_range(it->first);
		multimap<int,double>::iterator Charge_it = it_intern.first;
		if (Charge_it->second > charge_max)
		{
			charge_max = Charge_it->second;
			strip_max = Charge_it->first;
		}
	}
	return strip_max;
}

int clustr::MaxSoN()
{
	double SoN_max = 0;
	int strip_max;
	pair<multimap<int,double>::iterator, multimap<int,double>::iterator> it_intern; //New multimap iterator which is the return value of equal_range
	for (multimap<int,double>::iterator it = ElementsNoise.begin(); it != ElementsNoise.end(); it++) //loop over all elements of the cluster ElementsNoise multimap
	{
		it_intern=ElementsNoise.equal_range(it->first); //hands back the multimap iterator of the element (*it).first which is the strip number
		multimap<int, double>::iterator Charge_it = it_intern.first; //first element of this iterator stores the charge value
		multimap<int, double>::iterator Noise_it = next(Charge_it); //the next element (second and last) stores the noise value
		
		double SoN = Charge_it->second/Noise_it->second;
		if (SoN > SoN_max)
		{
			SoN_max = SoN;
			strip_max = Charge_it->first;
		}
	}
	//cout << "Maximum signal over noise is at strip " << strip_max << " with Signal/Noise of " << SoN_max;
	return strip_max;
}


clustr::~clustr() { }
