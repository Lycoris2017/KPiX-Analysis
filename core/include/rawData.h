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
	
	class Cycle{
	public:
		// Global channel num = kpix_index*2014+channel;
		Cycle(KpixEvent &event,
		      uint buc_level,
		      uint begin_ch,
		      uint end_ch,
		      bool isold);
		~Cycle(){};
		uint m_cyclenumber;
		uint16_t m_ts;
		uint m_b_level;
		
		// Do your own boosted map
		vector<uint> m_v_hashkeys_b[4];
		vector<uint> hashkeys(uint i){
			if (i<4)return m_v_hashkeys_b[i];
			else {
				printf("hashkeys error\n");
				exit(EXIT_FAILURE);
			}
		}

		bool m_has_adc;
		// Please make sure both ADC and TS are 16bit (see KpixSample)
		vector<uint16_t> m_v_adc_b[4];
		vector<uint16_t> vadc(uint i){
			if (i<4)return m_v_adc_b[i];
			else {
				printf("vadc error\n");
				exit(EXIT_FAILURE);
			}
		}
		vector<uint16_t> m_v_ts_b[4];
		
		bool m_has_fc;
		vector<double> m_v_fc_b[4];
		
		void RemovePed(std::vector<uint16_t> &ped_adc);
		void RemoveCM_CalFc(std::vector<uint16_t> &slopes, bool remove_adc); //-> remove the m_v_adc_*
		
		static void AddAdcBuf(Cycle&);
		static void ResetAdcBuf(){
			for(uint i; i<4;i++)s_buf_adc[i].clear(); };
		static void ResetPedDB();
		static void CalPed();
		//	private:
		static vector<uint> s_hashkeys[4];
		static vector<uint16_t> s_ped_adc[4];
		static unordered_map<uint, vector<uint16_t>> s_buf_adc[4];
	public:
		static void AddFcBuf(Cycle&);
		static void ResetFcBuf();
		static void ResetFcDB();
		static void CalNoise();
	private:
		static vector<double> s_noise_fc;
		static vector<vector<double>> s_buf_fc;
	public:
		static void loadCalib(const std::string&); // Done
		static unordered_map<uint, double> s_slopes_b0; //Done
		
		// Hash table
		static uint hashCode(uint kpix, uint channel){ return kpix*1024+channel;}
		static uint getKpix(uint hashcode){ return hashcode/1024;}
		static uint getChannel(uint hashcode){ return hashcode%1024;}

	};
	
	class rawData{
	public:
		rawData(){};
		~rawData(){};

		void loadFile(const std::string&);

	private:
		std::vector<Cycle> m_v_cycles;
	};


    

}
#endif
