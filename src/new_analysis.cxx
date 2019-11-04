//-----------------------------------------------------------------------------
// Mengqing Wu <mengqing.wu@desy.de>
// @ 2019-10-16
// Target: analysis template using new lycoris database classes
//-----------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string.h>

#include "rawData.h"
#include "DataRead.h"

using namespace std;


// Process the data
int main ( int argc, char **argv ) {
  if (argc == 1){
    printf("[Usage] ./new_analysis [input.dat]\n");
    return 0;
  }

  printf("[Info] You choose file %s\n", argv[1]);
  Lycoris::rawData db;
  db.setMaxCycles(100);
  db.setNBuckets(1);
  
  db.loadFile(argv[1]);
  //  cout<< "[dev] How many cycles? "  << db.GetNEvents() << std::endl;
  db.doRmPed();
  //  cout<< "debug, size of ped_adc: "<< Lycoris::Cycle::s_ped_adc[0].size()<< endl;
  Lycoris::Cycle::loadCalib("/home/lycoris-dev/workspace/kpix-analysis/data/calib_HG_201907T24.csv");
  //  db.doRmCMnoise();
  
  return 1;
 
}
