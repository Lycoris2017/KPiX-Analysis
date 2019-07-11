//-----------------------------------------------------------------------------
// File          : analysis.h
// Author        : Mengqing Wu (orig. Uwe, Kraemer & Ryan,Herbst) <mengqing.wu@desy.de>
// Created       : 08/16/2018
// Project       : KPiX Analysis
//-----------------------------------------------------------------------------
// Description :
// General analysis of KPiX data.
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Modification history :
// 16/08/2018: refactorize kpix analysis code (orig. Ryan Herbst, then modified by Uwe Kraemer)
//-----------------------------------------------------------------------------

// C++ std header
#include <stdexcept>
#include <vector>
// kpix header
#include "DataRead.h"
#include "KpixEvent.h"
#include "KpixSample.h"
#include "Data.h"

// extern header
//#include "KpixRead.h"

// root header
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"

// Class definition
namespace lycoris
{

  class ntupleMaker {

  protected:
    DataRead               _dataRead;  //kpix event classes used for analysis of binary date
    KpixEvent              _cycleevent;    //
    KpixSample             *_cyclesample;   //
    KpixSample::SampleType  _type;

    uint                   _tstamp;
    uint                   _kpix;
    uint                   _channel;
    uint                   _bucket;
    uint                   _value;
    uint                   _range;

    bool                   _debug;
    
  public:
    ntupleMaker(const char* binFile);
    ~ntupleMaker(){};
    void setDebug(bool debug){_debug = debug;};
    void CreateTree();
    void makeTreeEx();
    void loopKpix();
    void addBranch(TTree &tree){};

    typedef struct {
      uint tstamp; // bunchClkCount based
      uint value; // unit of 1/8 of the bunchClkCount == based on AcqClk.
    }kEXTTS;

    typedef struct {
      uint tstamp; // bunchClkCount based
      uint kpix;
      uint channel;
      uint bucket;
      uint range;
      uint ADC;
    }kDATA;

    typedef struct {
      uint CycleNr;
      uint CycleSampleCount;
      uint CycleTimeStamp;
    }CYCLEN;
    
  };

  ntupleMaker::ntupleMaker(const char* binFile){
    // Load the binary file and check throw exception if not exist
    _debug = false;
    try {
      bool atest = _dataRead.open(binFile);
      if (!atest){
	throw std::invalid_argument( std::string("Error opening binary data file") );
      }
      else
	std::cout << "Opened data file: " << binFile << std::endl;
    }
    catch ( const std::invalid_argument& e){
      std::cerr << "\nexception: " << e.what() << std::endl; 
    }
    
  }

  void ntupleMaker::loopKpix(){
    /* 
     * Test function for reading kpix data out
     */
    
    int eventcount=0;
    while ( _dataRead.next(&_cycleevent) ){
      eventcount++;

    }// Loop over kpix event
    
    _dataRead.close(); // TODO: other func after this func will have trouble!
    
    cout<< "In total, we have #" << eventcount << " events :)\n" << endl;
    return;
  }
  
  /* void ntupleMaker::addBranch(TTree &tree) { */

  /* } */
    
  void ntupleMaker::CreateTree(){

    TFile tfile("kpixTree.root","RECREATE","[dev] kpix ROOT file with histograms & trees");
    // Create a ROOT tree
    TTree *t1 = new TTree("Cycles","temperature etc /cycle");
    TTree *t3 = nullptr;
    TTree *t2 = nullptr;
    TH1F  *htest = nullptr;
    if (_debug){
      htest   = new TH1F("htest","CycleNr",7000,83000,90000);
      t2 = new TTree("HitEvent","Hit/range/bucket/channel/kpix, with cycle index");
      t3 = new TTree("ExtTstamp","external time stamp / cycle");
    }

    // Define some templates for the tree branch
    
    //typedef struct {Float_t x,y,z;} POINT;
 
    //std::map<std::string, std::vector<double>> dictTest;
    std::vector<kEXTTS> vtstamp;
    std::vector<kDATA>  vHits;
    
    static kDATA Hits;
    static kEXTTS exttstamp;
    static CYCLEN cyclen;
    //static POINT hit;

    std::vector<uint> kpixOnRun;
    int cyclecount=0;
    int nHits=0; // nData per cycle (Data is /bucket/channel/kpix)
    int nTrig=0; // nTrig per cycle
    int ev=0;
    uint CycleNr=0;
    // Create TBranches
    t1->Branch("tstamp",&vtstamp);
    t1->Branch("Hits",&vHits);
    //t1->Branch("dictTest",&dictTest);
   
    t1->Branch("cyclen", &cyclen, "CycleNr/i:CycleTimeStamp/i:CycleSampleCount/i");
    t1->Branch("kpixOnRun", &kpixOnRun);
    t1->Branch("cyclecount", &cyclecount, "cyclecount/I");
    t1->Branch("nHits", &nHits, "nHits/I");
    t1->Branch("nTrig", &nTrig, "nTrig/i");

    if (t2!=nullptr){
      t2->Branch("ev", &ev, "ev/I");
      t2->Branch("CycleNr", &CycleNr, "CycleNr/i");
      t2->Branch("Hits", &Hits, "tstamp/i:kpix/i:channel/i:bucket/i:range/i:ADC/i");
    }
    if (t3!=nullptr){
      t3->Branch("CycleNr", &CycleNr, "CycleNr/i");
      t3->Branch("exttstamp", &exttstamp, "tstamp/i:value/i");
    }
    /* // Here we loop over kpix event to fill the tree: */
    while ( _dataRead.next(&_cycleevent) ){
      kpixOnRun.clear();
      vtstamp.clear();
      vHits.clear();
      //dictTest.clear();
      nHits=0;
      nTrig=0;
      cyclecount++;
      
      if (htest!=nullptr)   htest->Fill( _cycleevent.eventNumber() );
      
      cyclen.CycleNr = _cycleevent.eventNumber();
      cyclen.CycleSampleCount = _cycleevent.count();
      cyclen.CycleTimeStamp = _cycleevent.timestamp();

      CycleNr = _cycleevent.eventNumber(); // t2
	
      for (uint x1=0; x1 < _cycleevent.count(); x1++ ){
    	// Get kpixSample == tracker event
    	_cyclesample = _cycleevent.sample(x1);
	
    	if ( _cyclesample->getEmpty() ){
    	  cout<<" [info] empty sample, jump over!"<<endl;
    	  continue;
    	}

	_type    = _cyclesample -> getSampleType();
	_tstamp  = _cyclesample -> getSampleTime();
	_value   = _cyclesample -> getSampleValue();// can be ADC(Data), TS(TimeStamp), etc
	_range   = _cyclesample -> getSampleRange();// for low/normal/high Gain dynamic change registered in KPiX
	_kpix    = _cyclesample -> getKpixAddress();
    	_channel = _cyclesample -> getKpixChannel();
    	_bucket  = _cyclesample -> getKpixBucket();
					
    	//if (!_kpix) cout << "debug: kpix address = " << _kpix <<endl;

	if (_type == KpixSample::Timestamp) {
	  nTrig++;
	  exttstamp.tstamp = _tstamp;
	  //exttstamp.value  = _value;
	  exttstamp.value = _cyclesample -> getBunchCount();
	  vtstamp.push_back(exttstamp);

	  //std::vector<double> _test;
	  //_test.push_back(1.9526);
	  //dictTest.insert(std::make_pair("123",_test));
	  if (t3!=nullptr)   t3->Fill();
	}
	
    	if ( _type == KpixSample::Data ){
    	  if ( std::find(kpixOnRun.begin(), kpixOnRun.end(), _kpix) == kpixOnRun.end() ) 
	    kpixOnRun.push_back(_kpix);
	  nHits++;
	  ev++;

	  Hits.tstamp  = _tstamp;
	  Hits.kpix    = _kpix;
	  Hits.channel = _channel;
	  Hits.bucket  = _bucket;
	  Hits.ADC     = _value;
	  Hits.range   = _range;
	  if (t2!=nullptr) t2->Fill();

	  vHits.push_back(Hits);
	  
    	}
      }

      t1->Fill();
     
    }// Loop over kpix event

    t1->Print();
    //t1->Write();
    tfile.Write();
    tfile.Close();
      
    _dataRead.close(); // TODO: other func after this func will have trouble!
    
    cout<< "In total, we have #" << cyclecount << " events :)\n" << endl;
        
  }


  //---------------------
  
  void ntupleMaker::makeTreeEx(){
    // do sth
    // Create a new ROOT binary machine independent file.
    // Note that this file may contain any kind of ROOT objects, histograms,trees
    // pictures, graphics objects, detector geometries, tracks, events, etc..
    // This file is now becoming the current directory.
    TFile hfile("htree.root","RECREATE","Demo ROOT file with histograms & trees");
    // Create some histograms and a profile histogram

    TH1F *hpx   = new TH1F("hpx","This is the px distribution",100,-4,4);
    TH2F *hpxpy = new TH2F("hpxpy","py ps px",40,-4,4,40,-4,4);
    TProfile *hprof = new TProfile("hprof","Profile of pz versus px",100,-4,4,0,20);
    // Define some simple structures
    typedef struct {Float_t x,y,z;} POINT;
    typedef struct {
       Int_t ntrack,nseg,nvertex;
       UInt_t flag;
       Float_t temperature;
    } EVENTN;
    static POINT point;
    static EVENTN eventn;
    // Create a ROOT Tree
    TTree *tree = new TTree("T","An example of ROOT tree with a few branches");
    tree->Branch("point",&point,"x:y:z");
    tree->Branch("eventn",&eventn,"ntrack/I:nseg:nvertex:flag/i:temperature/F");
    tree->Branch("hpx","TH1F",&hpx,128000,0);
    Float_t px,py,pz;
    static Float_t p[3];
    // Here we start a loop on 1000 events
    for ( Int_t i=0; i<1000; i++) {
       gRandom->Rannor(px,py);
       pz = px*px + py*py;
       Float_t random = gRandom->Rndm(1);
       // Fill histograms
       hpx->Fill(px);
       hpxpy->Fill(px,py,1);
       hprof->Fill(px,pz,1);
       // Fill structures
       p[0] = px;
       p[1] = py;
       p[2] = pz;
       point.x = 10*(random-1);;
       point.y = 5*random;
       point.z = 20*random;
       eventn.ntrack  = Int_t(100*random);
       eventn.nseg    = Int_t(2*eventn.ntrack);
       eventn.nvertex = 1;
       eventn.flag    = Int_t(random+0.5);
       eventn.temperature = 20+random;
       // Fill the tree. For each event, save the 2 structures and 3 objects
       // In this simple example, the objects hpx, hprof and hpxpy are slightly
       // different from event to event. We expect a big compression factor!
       tree->Fill();
    }
    // End of the loop
    tree->Print();
    // Save all objects in this file
    hfile.Write();
    // Close the file. Note that this is automatically done when you leave
    // the application.
    hfile.Close();

  }



  
}// namespace lycoris end
