#include "rawData.h"
#include <KpixEvent.h>
#include <KpixSample.h>
#include "Data.h"
#include "DataRead.h"

using namespace std;



Lycoris::rawData::rawData(){
  printf("rawData object constructor.\n");
}

Lycoris::rawData::~rawData() { }

void Lycoris::rawData::loadFile(const std::string& fname, bool isold){
  DataRead               dataRead;
  //  off_t                  fileSize;  //
  //  off_t                  filePos;   //
  KpixEvent              event;    //
  KpixSample             *sample;   //

  if ( ! dataRead.open(fname.c_str())  ) {
    cout << "Error opening data file " << fname << endl;
    return;
  }
  else printf("Open file %s\n", fname.c_str());
  cout << "[debug] " << dataRead.debug_ << endl;;
  while ( dataRead.next(&event) ){
    uint eventnumber, runtime;
    eventnumber = event.eventNumber();
    if (isold) runtime = event.timestamp();
    else       runtime = event.runtime();

    std::vector<Lycoris::rawTrack> tracks(4, rawTrack());

    for (uint ev=0; ev< event.count(); ev++){

      sample = event.sample(ev);

      uint kpix, channel, bucket, value;
      uint bunchClk, subCount, tstamp;
      
      if (sample->getSampleType() == KpixSample::Data){
	kpix    = sample->getKpixAddress();
	channel = sample->getKpixChannel();
	bucket  = sample->getKpixBucket();
	value   = sample->getSampleValue();
	tstamp  = sample->getSampleTime();
	Lycoris::rawHit hit(kpix, channel, bucket, value, tstamp);
	tracks.at(bucket).AddHit(hit);
	
      }
      
      if (sample->getSampleType() == KpixSample::Timestamp){
	tstamp   = sample->getSampleRuntime64(runtime);
	subCount = sample->getSubCount();
	bunchClk = sample->getBunchCount();

	double time = bunchClk + double(0.125 * subCount);
	
	std::cout<<"[dev] Trigger ts = " << time 
		 <<", runtime[ns] = " <<  tstamp*5 << std::endl;
      }
      

    }// finish loop over all data of one cycle
    
    m_rawdb[eventnumber] = tracks;

  }
  
  
}
