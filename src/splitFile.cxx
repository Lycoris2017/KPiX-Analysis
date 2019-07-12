//-----------------------------------------------------------------------------
// File          : printDat.cxx
// Author        : Mengqing Wu (origin from Ryan@SLAC)
// Created       : 01/10/2018
// Project       : simple evt loop to print data.
//-----------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <stdarg.h>
#include <math.h>
#include <fstream>
#include <string.h>

#include "TH1F.h"
#include "TFile.h"

#include "DataRead.h"
#include "KpixEvent.h"
#include "KpixSample.h"
#include "Data.h"

using namespace std;


// Process the data
int main ( int argc, char **argv ) {
  DataRead               dataRead;
  //off_t                  fileSize;
  //off_t                  filePos;
  KpixEvent              event;
  KpixSample             *sample;

  uint64_t               frameruntime;
  uint64_t               sampleruntime;

  KpixSample::SampleType type;

  uint                   kpix;
  uint                   channel;
  uint                   bucket;
  uint                   kpixeventcount;
  stringstream           tmp;
  ofstream					ofile("test.bin", ios::binary);
  
  
  // Data file is the first and only arg
  if ( argc != 2  ) {
    cout << "\nUsage: ./count data_file\n";
    return(1);
  }
  
  // Open data file
  if ( ! dataRead.open(argv[1])  ) {
    cout << "Error opening data file " << argv[1] << endl;
    return(1);
  }

  cout << "Opened data file: " << argv[1] << endl;
  //fileSize = dataRead.size();
  //filePos  = dataRead.pos();
  kpixeventcount = 0;
  cout << "\rReading File: 0 %" << flush;
  // Open root file
  tmp.str("");
  tmp<<"printDat.root";
  string outRoot = tmp.str();
  

  // read over kpix events again to fill the histograms:
  //dataRead.open(argv[1]);
  while ( dataRead.next(&event) ){
    cout << "DEBUG 1 : " << event.count() << endl;
    ofile.write((char*) &event, sizeof(KpixEvent));
  }
  dataRead.close(); // close file as we have looped through it and are now at the end

  cout << "Writing root plots to " << outRoot << endl;
  
  return 1;
	
}
