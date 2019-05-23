


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
	if (sensor == 0 || sensor == 5 || sensor == 4) // kpix side showing towards beam movement beam side  KPIX >  < Beam
	{
		if (sensor == 5) return 50*(-strip+1840); // sensor 2 and sensor 5 have no stereo angle
		else return 49.96954*(-strip+1840); // sensor 0,1, 3 and 4 have a 2 degree stereo angle means the y movement is slightly lower (50*cos(+-2/57.295))

	}
	else  // kpix side in direction of beam movement KPIX < < BEAM
	{
		if (sensor == 2) return 50*strip; // sensor 2 and sensor 5 have no stereo angle
		else return 49.96954*strip; // sensor 0,1, 3 and 4 have a 2 degree stereo angle means the y movement is slightly lower (50*cos(+-2/57.295))
	}
}

double xParameterSensor(double strip, int sensor)
{
	if (sensor == 0 || sensor == 5 || sensor == 4) // kpix side showing towards beam movement beam side  KPIX >  < Beam
	{
		if (sensor == 5) return 0*(-strip+920); // sensor 2 and sensor 5 have no stereo angle
		else return -1.744998*(-strip+920); // sensor 0 and 4 have a -2 degree stereo angle means the y movement is slightly lower (50*sin(-2/57.295))

	}
	else  // kpix side in direction of beam movement KPIX < < BEAM
	{
		if (sensor == 2) return 0*(strip-920); // sensor 2 and sensor 5 have no stereo angle
		else return 1.744998*(strip-920); // sensor 1 and 3 have a 2 degree stereo angle means the y movement is slightly lower (50*sin(2/57.295))
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
