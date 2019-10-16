#include "rawData.h"
#include <KpixEvent.h>
#include <KpixSample.h>
#include <Data.h>
#include <DataRead.h>

using namespace std;

rawData::rawData()
{
	kpix     = 9999;
	channel  = 9999;
	bucket   = 9999;
	range    = 9999;
	value    = 9999;
	time     = 9999;
	event    = 9999;
}

rawData::~rawData() { }

void rawData::loadFile(const std::string& fname, bool isold){
  DataRead               dataRead;
  off_t                  fileSize;  //
  off_t                  filePos;   //
  KpixEvent              event;    //
  KpixSample             *sample;   //

  if ( ! dataRead.open(fname.c_str())  ) {
    cout << "Error opening data file " << fname << endl;
    return;
  }
  while ( dataRead.next(&event) ){
    uint eventnumber, runtime;
    eventnumber = event.eventNumber();
    if (isold) runtime = event.timestamp();
    else       runtime = event.runtime();
    
    for (uint ev=0; ev< event.count(); ev++){
      sample = event.sample(ev);

      uint kpix, channel, bucket, value;
      uint bunchClk, subCount, timestamp;
      if (sample->getSampleType() == KpixSample::Data){
	kpix    = sample->getKpixAddress();
	channel = sample->getKpixChannel();
	bucket  = sample->getKpixBucket();
	value   = sample->getSampleValue();
	
	cout<<"[dev] kpix = "<<kpix<<", channel = " <<channel <<", bucket = " <<bucket <<"\n";
      }
      if (sample->getSampleType() == KpixSample::Timestamp){
	double time = bunchClk + double(0.125 * subCount);
	timestamp = sample->getSampleRuntime64(runtime);

	std::cout<<"[dev] Trigger ts = " << timestamp << std::endl;
      }
      

    }


  }
  
  
}
