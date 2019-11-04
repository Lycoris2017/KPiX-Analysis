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
double smallest_time_diff(std::vector<double>, int);

template <typename T>
T median(std::vector<T>* vec){
  /* Note by Mengqing: this func has already an optimal in memory/speed*/
  if (vec  == nullptr ){
      //cout << "Found a nullpointer" << endl;
    return 0;
  }
  else{
    if (vec->empty())  return 0;
    else{
      size_t n = vec->size() / 2;
      if (vec->size()%2 == 0){
	std::nth_element(vec->begin(), vec->begin()+n, vec->end());
	std::nth_element(vec->begin(), vec->begin()+n-1, vec->end());
	return (vec->at(n)+vec->at(n-1))/2;
      }
      else{
	std::nth_element(vec->begin(), vec->begin()+n, vec->end());
	return vec->at(n);
      }
    }
  }
}

template <typename T>
T MAD(std::vector<T>* v){
  if (v  == nullptr )
    {
      //cout << "Found a nullpointer" << endl;
      return 0;
    }
  else
    {
      if (v->empty())
	  return 0;
      else{
	T med = median(v);
	std::vector<T>* deviation;
	deviation = new std::vector<T>;
	
	for (auto const i:(*v)){
	  deviation->push_back(fabs(i - med));
	}
	return median(deviation);
      }
    }
}
#endif
