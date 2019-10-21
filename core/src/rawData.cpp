#include "rawData.h"
#include <KpixEvent.h>
#include <KpixSample.h>
#include "Data.h"
#include "DataRead.h"
#include <fstream>
#include <cstring>
#include <chrono>
#include <sstream>
using namespace std;


Lycoris::rawData::rawData(){
  printf("rawData object constructor.\n");

  //- open default geo
  string geo = getenv("KPIX_ANA");
  geo += "/data/plane_Geo_default.txt";
  loadGeo(geo);
string calib = getenv("KPIX_ANA");
  calib+= "/data/calib_HG_201907T24.csv";
  loadCalib(calib);
}

Lycoris::rawData::~rawData() { }

void Lycoris::rawData::loadCalib(const std::string& calib){
  // based on a csv file 

  printf(" loadCalib...\n");
  std::ifstream file (calib);
  if ( !file ) return;
  printf(" Open Calib file : %s\n", calib.c_str());
  
  std::string line;
  char delim=',';

  //m_calib_map.clear();
  while( std::getline(file, line) ){

    if (line.empty() ) continue;

    while (line[0]==' ')
      line.erase(0,1);

    if (!isdigit(line[0])) continue;

    std::stringstream is(line);
    std::string value;
    std::vector<std::string> vec;
    while (std::getline(is, value, delim))
      vec.push_back(value);
    if (vec.size()<4){
      printf("[ERROR] Missing value at");
      for (const auto &a : vec)
	cout << a << ' ';
      printf("\n");
    }
    else{
      auto kpix    = std::atoi(vec[0].c_str());
      auto channel = std::atoi(vec[1].c_str());
      auto bucket  = std::atoi(vec[2].c_str());
      float slope  = std::atof(vec[3].c_str());
      if (bucket==0)
	m_calib_b0.emplace(ChanKey(kpix,channel), slope);
    }
  }

  cout << " test: "<< m_calib_b0.at(ChanKey(11, 1000)) << endl;
  file.close();
  
}

void Lycoris::rawData::loadGeo(const std::string& geo){
  printf(" loadGeo...\n");
 
  ifstream file(geo);
  if (!file )  return;
  printf(" Open Geo file : %s\n",geo.c_str());

  std::string line;
  char delim=',';

  m_kpix2plane.clear();
  while ( std::getline(file, line) ){
    if (!line.empty())
      while (line[0] == ' ') line.erase(0,1);
        
    if (!line.empty() && isdigit(line[0])){
      std::stringstream is(line);

      std::string value;
      std::vector<uint> pair;
      while (std::getline(is, value, delim) )
	pair.push_back(atoi(value.c_str()));
      
      m_kpix2plane.insert(std::make_pair(pair.at(0),pair.at(1)));
    }
  }
  file.close();

}

void Lycoris::rawData::loadFile(const std::string& fname, bool isold){
  DataRead               dataRead;
  off_t                  fileSize;  //
  off_t                  filePos;   //
  uint                   lastPct;
  uint                   currPct;
  KpixEvent              event;    //
  KpixSample             *sample;   //

  if ( ! dataRead.open(fname.c_str())  ) {
    cout << "Error opening data file " << fname << endl;
    return;
  }
  else printf("Open file %s\n", fname.c_str());
  fileSize = dataRead.size();
  filePos = dataRead.pos();
  // Init
  currPct          	= 0;
  lastPct          	= 100;
  
  cout << "[debug] " << dataRead.debug_ << endl;;
  auto t_start = std::chrono::high_resolution_clock::now();
  
  while ( dataRead.next(&event) ){
    uint eventnumber, runtime;
    eventnumber = event.eventNumber();
    if (isold) runtime = event.timestamp();
    else       runtime = event.runtime();

    std::vector<Lycoris::rawTrack> tracks(1, rawTrack());
    
    for (uint ev=0; ev< event.count(); ev++){

      sample = event.sample(ev);

      uint kpix, channel, bucket, value;
      uint bunchClk, subCount, tstamp;
      double fc;
      
      if (sample->getSampleType() == KpixSample::Data){
	kpix    = sample->getKpixAddress();
	channel = sample->getKpixChannel();
	bucket  = sample->getKpixBucket();
	value   = sample->getSampleValue();
	tstamp  = sample->getSampleTime();

	//-currently only look at bucket0
	if (bucket!=0) continue;
	if (m_calib_b0.at(ChanKey(kpix, channel)))
	  fc  =  value/m_calib_b0.at(ChanKey(kpix, channel));
	
	Lycoris::rawHit hit(kpix, channel, bucket, value, tstamp, fc);
	tracks.at(bucket).AddHit(kpix, hit);
	//tracks.at(bucket).AddHit(hit);

      }
      
      if (sample->getSampleType() == KpixSample::Timestamp){
	tstamp   = sample->getSampleRuntime64(runtime);
	subCount = sample->getSubCount();
	bunchClk = sample->getBunchCount();

	double time = bunchClk + double(0.125 * subCount);
	
	//std::cout<<"[dev] Trigger ts = " << time 
	//	 <<", runtime[ns] = " <<  tstamp*5 << std::endl;
      }
      

    }// finish loop over all data of one cycle
    
    m_rawdb[eventnumber] = tracks;

    filePos  = dataRead.pos();
    currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
    if ( currPct != lastPct ) {
      cout << "\rReading File for noise determination: " << currPct << " %      " << flush;
      lastPct = currPct;
    }
  }
  
  auto t_end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds> (t_end - t_start).count();
  cout << " Time[s]: " << duration
       << endl;

  dataRead.close();
}


void Lycoris::rawData::Pedestals() const{



}
