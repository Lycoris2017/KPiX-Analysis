#ifndef __RAWDATA_H__
#define __RAWDATA_H__
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <stdlib.h>
#include <iterator>
#include <map>
#include <vector>

/*
 Mengqing <mengqing.wu@desy.de>
 @ 2019-10-17
 Currently for external trigger only
 */
namespace Lycoris
{
  uint INIT=9999;
  // Loading a calication and translate ADC function can be added here also.
  
  class rawHit{
  public:
  rawHit():
    m_kpix(9999), m_channel(9999), m_bucket(9999), m_charge_adc(9999),m_ts(9999){}
    
  rawHit(uint kpix, uint channel, uint bucket, uint adc, uint ts ):
    m_kpix(kpix), m_channel(channel), m_bucket(bucket), m_charge_adc(adc),m_ts(ts){}
    
  private:
    uint m_kpix;
    uint m_channel;
    uint m_bucket;
    uint m_charge_adc;
    uint m_charge_fc;  // need calibration data
    uint m_ts;
    double m_pedestal; // over all cycle, need calculation
    double m_cm_noise; // after pedestal, need calculation
  };

  /* class localPlane{ */
  /* public: */
  /*   // require geometry input, say which two ASICs to which layer */
  /*   void loadGeo(const std::string&); */

  /*   // make a map to input for PacMan clustering algo */
  /*   // Strip, Noise/Charge */
  /*   std::map<int, double> makeChargeMap(); */
  /*   std::map<int, double> makeNoiseMap(); */
    
  /* private: */
  /*   typedef std::vector<rawHit> hits_t; */
  /*   hits_t m_hits; */
  /*   uint m_plane_id; */
    
  /* }; */
  
  
  class rawTrack{
  public:
    rawTrack(){};
    ~rawTrack(){};
    uint GetNPlanes();
    uint GetNKpixs();
    void AddHit(rawHit hit){m_hits.push_back(hit);};
    void MakePlanes(); // make plane out of all hits
    
  private:
    typedef std::vector<rawHit> hits_t;
    hits_t m_hits;

  };
    
  
  class rawData {
  public:
    // Constructor
    rawData();	
    // Deconstructor
    ~rawData();
    
    void loadFile(const std::string&, bool isold = false);
    size_t GetNEvents(){ return m_rawdb.size(); };
  private:
    typedef std::map<uint, std::vector<rawTrack> > rawDB_t;
    rawDB_t m_rawdb;
    
  } ;
  
}

#endif
