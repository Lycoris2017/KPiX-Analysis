#ifndef __TBFunctions_H__
#define __TBFunctions_H__
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <math.h> /* fabs */
#include <iostream>

using namespace std;

double yParameter(double, int);
double yParameterSensor(double, int);
double xParameterSensor(double, int);
double median(vector<double>*);
double MAD(vector<double>*);

#endif
