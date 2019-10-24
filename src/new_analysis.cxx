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
  db.loadFile(argv[1]);
  //  cout<< "[dev] How many cycles? "  << db.GetNEvents() << std::endl;
  
  return 1;
 
}
