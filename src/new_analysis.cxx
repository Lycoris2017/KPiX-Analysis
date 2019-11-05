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

#include "kpix_left_and_right.h"

#include "TH1F.h"
using namespace std;


// Process the data
int main ( int argc, char **argv ) {
  if (argc == 1){
    printf("[Usage] ./new_analysis [input.dat]\n");
    return 0;
  }

  printf("[Info] You choose file %s\n", argv[1]);
  Lycoris::rawData db;
  //db.setMaxCycles(1000);
  db.setNBuckets(1);
  db.loadFile(argv[1]);
  
  cout<< "[dev] How many cycles? "  << db.getNCycles() << std::endl;
  db.loadCalib("/home/lycoris-dev/workspace/kpix-analysis/data/calib_HG_201907T24.csv");

  db.doRmPedCM();
  uint b = 1;
  Lycoris::Cycle::CalNoise(b);

  
  TH1F *noise_distribution_sensor = new TH1F("h1", "noise_distribution; Noise(fC);   #channels", 100,-0.005, 0.995);

  
  return 1;
 
}
