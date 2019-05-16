


#include "TBFunctions.h"


using namespace std;


double yParameter(double strip, int kpix)
{

	if (kpix%2 == 0) // if left kpix
	{
		if (kpix == 6 || kpix == 2 || kpix == 4) //if lower kpix
		{
			return 50*strip;
		}
		else
		{
			return 50*(-strip+1840);
		}
	}
	else  // if right kpix
	{
		if (kpix == 1 || kpix == 9 || kpix == 11) // if lower kpix
		{
			return 50*(-strip+1840);
		}
		else
		{
			return 50*strip;
		}
		
	}
}

double yParameterSensor(double strip, int sensor)
{

	if (sensor == 0 || sensor == 3 || sensor == 4) // kpix side showing towards beam movement beam side  KPIX >  < Beam
	{
		return 50*(-strip+1840);
	}
	else  // kpix side in direction of beam movement KPIX < < BEAM
	{
		return 50*strip;
	}
}

double median(vector<double> &v)
{
    size_t n = v.size() / 2;
    if (n != 0)
    {
		if (v.size()%2 == 0)
		{
			nth_element(v.begin(), v.begin()+n, v.end());
			nth_element(v.begin(), v.begin()+n-1, v.end());
			return (v[n]+v[n-1])/2;
		}
		else
		{
			nth_element(v.begin(), v.begin()+n, v.end());
			return v[n];
		}
	}
	else return 0;
}

double MAD(vector<double> &v)
{
	double med = median(v);
	vector<double> deviation;
	for (auto const i:v)
	{
		deviation.push_back(fabs(i - med));
	}
	return median(deviation);
}
