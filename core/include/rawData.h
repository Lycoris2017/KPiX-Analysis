#ifndef __RAWDATA_H__
#define __RAWDATA_H__
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <stdlib.h>
#include <map>
#include <vector>
#include <unordered_map>

/*
 Mengqing <mengqing.wu@desy.de>
 @ 2019-10-17
 Currently for external trigger only
 */

namespace Lycoris
{
  const uint UINT_INIT=9999;
  
 
  // Loading a calication and translate ADC function can be added here also.
  
  class rawHit{
  public:
  rawHit():
    m_kpix(9999), m_channel(9999), m_bucket(9999), m_charge_adc(9999),m_ts(9999), m_charge_fc(0.0){}
    
  rawHit(uint kpix, uint channel, uint bucket, uint adc, uint ts, double fc):
    m_kpix(kpix), m_channel(channel), m_bucket(bucket), m_charge_adc(adc),m_ts(ts), m_charge_fc(fc){}
    
  private:
    uint m_kpix;
    uint m_channel;
    uint m_bucket;
    uint m_charge_adc;
    uint m_ts;
    double m_charge_fc;  // need calibration data
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
    /* Track is defined as per bucket, per event/cycle */
  public:
    rawTrack(){};
    ~rawTrack(){};
    //uint GetNPlanes(){return m_planes.size() };
    uint GetNKpixs() {return m_kpixs.size(); };
    void AddHit(rawHit hit){m_hits.push_back(hit);};
    void AddHit(uint kpix, rawHit hit){
      if (m_kpixs.count(kpix)) {// found key
	m_kpixs.at(kpix).push_back(hit);
      }
      else{
	hits_t _hits;
	_hits.push_back(hit);
	m_kpixs.emplace(kpix, _hits);
      }
      
    };

    /* class Key{ */
    /*   uint m_plane_id; */
    /*   uint m_kpix_id; */
    /* public: */
    /* Key(): m_plane_id(9999), m_kpix_id(9999) {} */
    /* Key(uint plane_id, uint kpix_id): m_plane_id(plane_id), m_kpix_id(kpix_id){} */
    /* } */
    
  private:
    typedef std::vector<rawHit> hits_t;
    hits_t m_hits;
    typedef std::unordered_map<uint, hits_t> hitCollections_t;
    hitCollections_t m_kpixs;
    //typedef std::unordered_map<Key, hits_t> planes_t;
    //planes_t m_planes;
    
  };
    
  
  class rawData {
  public:
    // Constructor
    rawData();	
    // Deconstructor
    ~rawData();
    
    void loadFile(const std::string&, bool isold = false);
    void loadGeo(const std::string&);
    void loadCalib(const std::string&);
    
    size_t GetNEvents() const{ return m_rawdb.size(); };
    uint GetPlaneIndex(uint kpix) const {return m_kpix2plane.at(kpix); };
    void Pedestals() const;
  private:
    // think about if rawDB_t can be changed to unordered_map or not to make sw faster.
    
    typedef std::unordered_map<uint, std::vector<rawTrack> > rawDB_t;
    rawDB_t m_rawdb;
    std::unordered_map<uint, uint> m_kpix2plane;

  public:
    class ChanKey{
    public:
      uint kpix;
      uint channel;
    ChanKey():kpix(9999), channel(9999) {}
    ChanKey(uint _kpix, uint _channel): kpix(_kpix), channel(_channel) {}
      
      bool operator ==( const ChanKey &k ) const {
	return kpix == k.kpix && channel == k.channel;
      }
      bool operator <(const ChanKey &k) const{
      	return kpix < k.kpix || ( kpix == k.kpix && channel < k.channel );
      }
    };
    
    
  private:
  typedef std::map<ChanKey, float> calib_map_t;
    calib_map_t m_calib_b0;
      
  } ;
  
}

#endif
