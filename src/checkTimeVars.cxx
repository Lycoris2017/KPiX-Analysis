//-----------------------------------------------------------------------------
// Mengqing Wu <mengqing.wu@desy.de>
// @ 2019-10-02
//
// Target: simple code to check different timestamps from Lycoris DAQ
//-----------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <stdarg.h>
#include <KpixEvent.h>
#include <KpixSample.h>
#include <Data.h>
#include <DataRead.h>
#include <math.h>
#include <fstream>
#include <string.h>

using namespace std;

int main( int argc, char **argv){
	stringstream   tmp;

	DataRead dataRead;
	if (! dataRead.open(argv[1])) {
		cout << "Error open data file "<< argv[1] << endl;
		return(1);
	}
	FILE* outcsv;
	tmp.str("");
	tmp << argv[1] << ".timing.csv" ;
	outcsv = fopen(tmp.str().c_str(),"w");
	printf("Write timing to file: %s\n", tmp.str().c_str());
	cout << "Open data file: "<< argv[1] << endl;
	fprintf(outcsv, "eventnumber, framruntime, frameruntime_ns, sampleruntime, sampleruntime_ns\n");

	cout << "Open data file: "<< argv[1] << endl;
	cout << "\rReading File: 0 %" << flush;
	// Init
	off_t fileSize = dataRead.size();
	off_t filePos  = dataRead.pos();
	uint currPct = 0;
	uint lastPct = 100;

	KpixEvent event; // cycle
	KpixSample *sample;  // min unit of data structure
	uint evcount=0;
	uint frameruntime=0;
	uint sampleruntime=0;
	uint eventnumber = 0;
	while (dataRead.next(&event)){
		evcount++;
		frameruntime = event.runtime(); //to print
		eventnumber  = event.eventNumber(); // to print
		for (uint ss=0; ss< event.count(); ss++){
			sample = event.sample(ss);
			if (sample->getEmpty()) continue;

			if ( (sample->getSampleType()) == KpixSample::Timestamp){
				sampleruntime = sample->getSampleRuntime64(frameruntime); // to print

				fprintf(outcsv, "%d,%d,%d,%d,%d\n", eventnumber, frameruntime, frameruntime*5, sampleruntime, sampleruntime*5);
			}
			
		}

		// show fancy data read bar
		filePos = dataRead.pos();
		currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
		if (currPct != lastPct){
			cout << "\rReading File: " << currPct << " %      " << flush;
			lastPct = currPct;
		}
	}
	cout << "\n"; // new line for flush
	
	return(0);
}
