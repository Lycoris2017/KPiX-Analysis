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
  bool                   doMiniAna=false; 
  bool                   printalot=false;
  bool                   print10evts=false;
  bool                   print10Sam=false;
  DataRead               dataRead;
  //off_t                  fileSize;
  //off_t                  filePos;
  KpixEvent              event;
  KpixSample             *sample;

  uint64_t               frameruntime;
  uint64_t               sampleruntime;
  
  bool                   bucketFound[32][1024][4];
  bool                   chanFound[32][1024];
  bool                   kpixFound[32];
  KpixSample::SampleType type;

  uint                   kpix;
  uint                   channel;
  uint                   bucket;
  uint                   kpixeventcount;
  double                 weight; // normalize to total kpix acq. cycles
    
  stringstream           tmp;
  TFile                  *rFile;
  TH1F*                  channel_entries[32][5]; //entries per bucket for each kpix, [kpix][bucket], bucket=5 is the bucket inclusive one;
  
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
  rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
  rFile->cd(); // move into root folder base
  
  //////////////////////////////////////////
  // Read Data - init
  /////////////////////////////////////////
  // check over how many kpixes w/ how many channels connected  
  printalot = false, print10evts=true, print10Sam=false;
  
  while ( dataRead.next(&event) ) {
    kpixeventcount++;

    /*if (kpixeventcount < 10 && print10evts) {
      cout << "[dev] 'KpixEvent' Level:" << endl;
      cout << "\tEventNumber = " << event.eventNumber()<<endl;
      cout << "\tTimeStamp = " << event.timestamp() <<endl;
      cout << "\tCount = " << event.count() <<endl;
      }*/
    //    cout << "Runtime : "<< event.runtime() << endl;
    frameruntime  = event.runtime();
    
    int DataSampleCount=0, TempSampleCount=0;
    int tenSamples=0;
    for (uint x1=0; x1 < event.count(); x1++) {
      //// Get sample
      sample  = event.sample(x1);
      if (sample->getEmpty()) {
	      //cout<<" [info] empty sample, jump over!"<<endl;
	      continue;
      }
      
      kpix    = sample->getKpixAddress();
      channel = sample->getKpixChannel();
      bucket  = sample->getKpixBucket();
      type    = sample->getSampleType();
      
      if ( type == KpixSample::Data ){
	      DataSampleCount++;
	      
	      kpixFound[kpix]          = true;
	      chanFound[kpix][channel] = true;
	      bucketFound[kpix][channel][bucket] = true;
	      
	      if ( print10Sam &&  tenSamples < 11 ){
		      //cout << kpixeventcount << endl;
		      cout<<"[dev] kpix = "<<kpix<<", channel = " <<channel <<", bucket = " <<bucket <<"\n";
		      tenSamples++;
	      }
      }
      
      if ( type == KpixSample::Temperature){
	      TempSampleCount++;
	      if ( print10Sam &&  tenSamples < 11 ){
		      cout<<"[*]   This is Temperature!" << endl;
		      tenSamples++;
	      }
      }
      if (type == KpixSample::Timestamp){
	if (frameruntime!=0){
	  sampleruntime = sample->getSampleRuntime64(frameruntime);
	}
	      else
		      cerr<< "Error: frameruntime is ZEROs!"<< endl;
      }
      
      kpixFound[0] = false; // in any case, kpix=0 is a virtual index
    }
    
    
    if (kpixeventcount < 10 && print10evts){
	    cout << "\tDataSampleCount = " << DataSampleCount <<endl;
	    cout << "\tTempSampleCount = " << TempSampleCount <<endl;
    }

  }
  
  dataRead.close();

  cout<< "In total, we have #"
      << dec << kpixeventcount
      << " events :)\n"<<endl;

  if (kpixeventcount>0) weight = 1.0/kpixeventcount;
  else {
    cout << "Error: kpix acq. cycles = "<< kpixeventcount<<endl;
    return 0;
  }

  //////////////////////////////////////////
  // Read Data
  /////////////////////////////////////////
  doMiniAna=false;
  if (!doMiniAna)  return 0;
  //  Initialize all histograms
  vector<uint> numkpix_vec;
  for (kpix=0; kpix<32; kpix++) {
    if (kpixFound[kpix]){
      numkpix_vec.push_back(kpix);

      tmp.str("");
      tmp << "Channel_entries_k" << kpix << "_allBucket";
      channel_entries[kpix][4] = new TH1F(tmp.str().c_str(), "", 1024, -0.5, 1023.5);

      for (bucket=0; bucket<4; bucket++){
	tmp.str("");
	tmp << "Channel_entries_k" << kpix << "_b" << bucket;
	channel_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel Entries; kpix channel addr.; #OfEvts/#acq.cycles", 1024, -0.5, 1023.5);
      }
    }
  }
  cout << "We have Number of Kpix ON = " << numkpix_vec.size()
       << "\n they are: ";
  for (auto& numkpix:numkpix_vec) cout<<numkpix<<" ";
  cout<<"\n\n";

  // read over kpix events again to fill the histograms:
  dataRead.open(argv[1]);
  while ( dataRead.next(&event) ){
    for (uint x2=0; x2<event.count(); x2++){
      //// Get sample
      sample  = event.sample(x2);  // check event subtructure
      if (sample->getEmpty()) continue; // if empty jump over
      
      kpix    = sample->getKpixAddress();
      channel = sample->getKpixChannel();
      bucket  = sample->getKpixBucket();
      //value   = sample->getSampleValue();
      type    = sample->getSampleType();
      //tstamp  = sample->getSampleTime();
      //range   = sample->getSampleRange();
      
      if (type == KpixSample::Data) {
	channel_entries[kpix][bucket]->Fill(channel, weight);
	channel_entries[kpix][4]->Fill(channel, weight);
      }
    }
  }
  dataRead.close(); // close file as we have looped through it and are now at the end

  cout << "Writing root plots to " << outRoot << endl;
  rFile->Write();
  rFile->Close();
  
  return 1;
	
}
