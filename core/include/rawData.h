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
#include <KpixEvent.h>

/*
 Mengqing <mengqing.wu@desy.de>
 @ 2019-10-17
 @ finalized 2019-10-24
 Currently for external trigger only
 */

namespace Lycoris{

  constexpr uint G_BUCKET_HASH_UNIT { 100000 };
  
	class Cycle{
	public:
		// Global channel num = kpix_index*2014+channel;
		Cycle(KpixEvent &event,
		      uint nbuckets = 1,
		      uint begin_ch = 0,
		      uint end_ch   = 1023,
		      bool isold    = false);
		~Cycle(){};
		uint m_cyclenumber;
		uint16_t m_ts;
		uint m_nbuckets;
		double m_cm_noise;
		
		// Do your own boosted map
		vector<uint> m_v_hashkeys_b[4];
		const vector<uint>& hashkeys(uint i) const{
			if (i<4)return m_v_hashkeys_b[i];
			else {
				printf("hashkeys error\n");
				exit(EXIT_FAILURE);
			}
		}

		bool m_has_adc;
		// Please make sure both ADC and TS are 16bit (see KpixSample)
		std::vector<uint16_t> m_v_adc_b[4];
		const vector<uint16_t>& vadc(uint b) const{
			if (b<4) return m_v_adc_b[b];
			else {
				printf("vadc error\n");
				exit(EXIT_FAILURE);
			}
		}
		vector<uint16_t> m_v_ts_b[4];
		
		bool m_has_fc;
		vector<double> m_v_fc_b[4];

		void RemovePed(std::unordered_map<uint, std::pair<uint, uint>> &, bool);
		void RemoveCM_CalFc(std::unordered_map<uint, double> &, bool ); //-> remove the m_v_adc_*
		
		static void AddAdcBuf(Cycle&);
		static void ResetAdcBuf(){
		  s_buf_adc.clear();
		  printf("[info] Static map ADC Buffer cleared.\n");
		}
		static void ResetPedMad(){
		  s_ped_adc_mad.clear();
		  printf("[info] Static map Ped ADC cleared");
		}
		static void CalPedMad(uint nbuckets =1);
		//	private:
		static unordered_map<uint, vector<uint16_t>> s_buf_adc;
		static unordered_map<uint, std::pair<uint, uint> > s_ped_adc_mad;

	public:
		static void AddFcBuf(Cycle&);
		static void ResetFcBuf();
		static void ResetNoise();
		static void CalNoise();
	private:
		static vector<double> s_noise_fc;
		static vector<vector<double>> s_buf_fc;
	public:
		static void loadCalib(const std::string&); // Done
		static unordered_map<uint, double> s_slopes_b0; //Done
		
		// Hash table
		static uint hashCode(uint kpix, uint channel){ return kpix*1024+channel;}
		static uint hashCode(uint kpix, uint channel, uint bucket){return bucket*G_BUCKET_HASH_UNIT + kpix*1024+channel;}
		static uint getKpix(uint hashcode){ return hashcode/1024;}
		static uint getChannel(uint hashcode){ return hashcode%1024;}
		static uint getBucketRm(uint hashcode, uint bucket) { return hashcode - bucket*G_BUCKET_HASH_UNIT;}

	};
	
	class rawData{
	public:
	rawData(): m_nbuckets(1), m_nmax(0){};
	  ~rawData(){};

	  void loadFile(const std::string&);
	  void doRmPed(bool rmMAD0=true); // TBD: MAD0 not yet 
	  void doRmCMnoise(bool rmAdc=true);
	  
	  void setNBuckets(uint level){
	    if (level<1 || level>4) return;
	    else m_nbuckets = level;
	  }
	  uint getNBuckets(){return m_nbuckets;}
	  void setMaxCycles(uint nmax){
	    printf("[info] Set max cycles to load: %d\n",nmax);
	    m_nmax = nmax;
	  }
	  
	private:
		std::vector<Cycle> m_v_cycles;
		uint m_nbuckets;
		uint m_nmax;
	};


    

}
#endif
