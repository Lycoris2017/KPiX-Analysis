#ifndef __TBFunctions_H__
#define __TBFunctions_H__
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <math.h> /* fabs */
#include <iostream>

//using namespace std;

double yParameter(double, int);
double yParameterSensor(double, int);
double xParameterSensor(double, int);
double MAD(std::vector<double>*);
double smallest_time_diff(std::vector<double>, int);

template <typename T>
T median(std::vector<T>* v){
  /* Note by Mengqing: this func has already an optimal in memory/speed*/
  if (v  == nullptr )
    {
      //cout << "Found a nullpointer" << endl;
      return 0;
    }
  else
    {
      if (v->empty())
	{
	  return 0;
	}
      else
	{
	  size_t n = v->size() / 2;
	  if (v->size()%2 == 0)
	    {
	      std::nth_element(v->begin(), v->begin()+n, v->end());
	      std::nth_element(v->begin(), v->begin()+n-1, v->end());
	      return (v->at(n)+v->at(n-1))/2;
	    }
	  else
	    {
	      std::nth_element(v->begin(), v->begin()+n, v->end());
	      return v->at(n);
	    }
	}
    }
}

#endif
