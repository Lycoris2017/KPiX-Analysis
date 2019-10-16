#ifndef __RAWDATA_H__
#define __RAWDATA_H__
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iterator>

using namespace std;

class rawData
{
 public:
  // Constructor
  rawData();	
  // Deconstructor
  ~rawData();
		
  int kpix;
  int channel;
  int bucket;
  int range;
  double value;
  double time;
  int event;

  void loadFile(const std::string&, bool isold = false);
	
} ;



#endif
