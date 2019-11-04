#include <KpixEvent.h>
#include <KpixSample.h>
#include <fstream>
#include <cstring>
#include <chrono>
#include <sstream>
#include <cassert>

#include "rawData.h"
#include "TBFunctions.h"
#include "Data.h"
#include "DataRead.h"

using namespace std;
using namespace Lycoris;

std::unordered_map<uint, double> Cycle::s_slopes_b0;
//std::vector<uint> Cycle::s_hashkeys[4];
std::vector<double> Cycle::s_noise_fc;
std::vector<vector<double>> Cycle::s_buf_fc;

std::unordered_map< uint, uint16_t> Cycle::s_ped_adc;
std::unordered_map< uint, vector<uint16_t>> Cycle::s_buf_adc;

Cycle::Cycle(KpixEvent &event, uint nbuckets,
                      uint begin_ch, uint end_ch,
                      bool isold ){

	// Declare member variables first
	m_nbuckets = nbuckets;
	m_has_adc = false;
	m_has_fc  = false;
	m_cyclenumber = event.eventNumber();
	if (isold) m_ts = event.timestamp();
    else       m_ts = event.runtime();
	

    KpixSample *sample;   //
    for (uint ev=0; ev< event.count(); ev++){

      sample = event.sample(ev);

      uint kpix, channel, bucket, value;
      uint bunchClk, subCount, tstamp;
      kpix    = sample->getKpixAddress();
      channel = sample->getKpixChannel();
      bucket  = sample->getKpixBucket();
      value   = sample->getSampleValue();
      tstamp  = sample->getSampleTime();

      if (sample->getSampleType() == KpixSample::Data){
	      if ( bucket >= m_nbuckets ) continue;
	      if ( channel < begin_ch || channel > end_ch)  continue;
	      auto key = hashCode(kpix, channel);

	      m_v_hashkeys_b[bucket].push_back(key);
	      m_v_adc_b[bucket].push_back(value);
	      m_v_ts_b[bucket].push_back(tstamp);
      }
      
      if (sample->getSampleType() == KpixSample::Timestamp){
	      double time;
	      if (isold)
		      time = tstamp + double(value * 0.125);
	      else{
		      tstamp   = sample->getSampleRuntime64(m_ts);
		      subCount = sample->getSubCount();
		      bunchClk = sample->getBunchCount();
		      time = bunchClk + double(0.125 * subCount);
	      }
	      // TBD: Do nothing now, add it to data structure later.
	      
	      //std::cout<<"[dev] Trigger ts = " << time 
	      //	 <<", runtime[ns] = " <<  tstamp*5 << std::endl;
      }
      

    }// finish loop over all data of one cycle
    for (uint i=0; i<4; i++){
      //printf("adc vec size : %d\n", m_v_adc_b[i].size() );
      m_has_adc = m_has_adc || m_v_adc_b[i].size();
	    
	    // **Exception Handler - Start
	    try{
		    if (m_v_hashkeys_b[i].size() != m_v_adc_b[i].size())
			    throw m_cyclenumber;
	    }
	    catch (uint e){
		    cout << "ERROR: Cycle Event " << e << " has channel_keys != its value. \n";
	    }
	    // **Exception Handler - End
    }
}


void Cycle::loadCalib(const std::string& fname){
  /* 
     Read Calib from a csv file
     cols needs to be:
     kpix >> channel >> bucket >> slope
   */
	s_slopes_b0.clear();
	printf(" loadCalib...\n");
	std::ifstream file (fname);
	if ( !file ) return;
	printf(" Open Calib file : %s\n", fname.c_str());
	
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
				Cycle::s_slopes_b0.emplace(hashCode(kpix,channel), slope);
		}
	}
	
	cout << " test: "<< Cycle::s_slopes_b0.at(hashCode(11, 1000)) << endl;
	file.close();
	
}

// void Lycoris::rawData::loadGeo(const std::string& geo){
//   printf(" loadGeo...\n");
 
//   ifstream file(geo);
//   if (!file )  return;
//   printf(" Open Geo file : %s\n",geo.c_str());

//   std::string line;
//   char delim=',';

//   m_kpix2plane.clear();
//   while ( std::getline(file, line) ){
//     if (!line.empty())
//       while (line[0] == ' ') line.erase(0,1);
        
//     if (!line.empty() && isdigit(line[0])){
//       std::stringstream is(line);

//       std::string value;
//       std::vector<uint> pair;
//       while (std::getline(is, value, delim) )
// 	pair.push_back(atoi(value.c_str()));
      
//       m_kpix2plane.insert(std::make_pair(pair.at(0),pair.at(1)));
//     }
//   }
//   file.close();

// }

void rawData::loadFile(const std::string& fname){
  DataRead               dataRead;
  off_t                  fileSize;  //
  off_t                  filePos;   //
  uint                   lastPct;
  uint                   currPct;
  KpixEvent              event;    //
  uint                   ncys=0;
  
  m_v_cycles.clear();

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
  
  //cout << "[debug] " << dataRead.debug_ << endl;;
  auto t_start = std::chrono::high_resolution_clock::now();

  Cycle::ResetAdcBuf();
  while ( dataRead.next(&event) ){
    if (m_nmax!=0 && ncys > m_nmax) break;
    ncys++;

    Cycle cy(event, m_nbuckets );
    // For pedestal calculation:
    Cycle::AddAdcBuf(cy); 
    // !!Once after std::move, you should never use this object
    m_v_cycles.push_back(std::move(cy));
    
    //---- fancy bar - Start 
    filePos  = dataRead.pos();
    currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
    if ( currPct != lastPct ) {
      cout << "\rLoading File for Cycle processing: " << currPct << " %      " << flush;
      lastPct = currPct;
    }
    //---- fancy bar - End 
    
  }
  
  auto t_end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds> (t_end - t_start).count();
  cout << " Time[s]: " << duration
       << endl;

  // //-- Debug print: Start
  // cout << "# of chan    :  " << Cycle::s_buf_adc[0].size() << endl;
  // cout << "# of cy in dt:  " << m_v_cycles.size() << endl;

  // uint count=0;
  // vector<uint> chan_cys;
  // for (const auto& x : Cycle::s_buf_adc[0]){
  // 	  if (x.second.size()!=m_v_cycles.size() ){
  // 		  count++;
  // 		  chan_cys.push_back(x.second.size());
  // 		  /*cout << "chan "
  // 		       << Cycle::getKpix(x.first) << "."
  // 		       << Cycle::getChannel(x.first) 
  // 		       << ": "
  // 		       << x.second.size() << endl;
  // 		  */
  // 	  }
  // }
  // sort( chan_cys.begin(), chan_cys.end() );
  // chan_cys.erase( unique( chan_cys.begin(), chan_cys.end() ), chan_cys.end() );
  // cout << "# of cy in chan: ";
  // for (const auto& x: chan_cys)
  // 	  cout << x << ", ";
  // cout << "\n";
  // cout << "# of such chan? " << count << endl;
  // //-- Debug print: End
  dataRead.close();
}

//- Fill s_ped_adc unordered_map<uint, uint16_t>
void Cycle::CalPed(uint nbuckets){
  
  if (s_buf_adc.empty()) return;
  for( auto &buf : s_buf_adc){
    
    // buf.first, buf.second -> median
    auto value = median(&buf.second);
    s_ped_adc.insert(std::make_pair(buf.first, value));
  }
  printf("Pedestal calculation finished\n");

  // after pedestal calculation, clean the buffer:
  Cycle::ResetAdcBuf();
}

//- Fill s_buf_adc unordered_map<uint, vec<uint16_t>>
void Cycle::AddAdcBuf(Cycle& cy){
  for (uint bb = 0; bb<cy.m_nbuckets; bb++){
    auto target = &s_buf_adc;
    auto kk = cy.hashkeys(bb);
    auto vv = cy.vadc(bb);
    
    /*if (vv.size() == 0 )
      printf("Cycle is empty: %d\n", cy.m_cyclenumber);
    */
    
    if(target->empty()){
      printf("Empty for evt: %d \n", cy.m_cyclenumber);
      // first event, init the vectors
      assert(kk.size() == vv.size());
      
      for( size_t cc = 0; cc < vv.size(); ++cc){
	std::vector<uint16_t> vec;
	vec.push_back(vv.at(cc));
	// key is kk.at(cc), value is vv.at(cc)
	auto key = Cycle::hashCode(getKpix(kk.at(cc)),
				   getChannel(kk.at(cc)),
				   bb);
	target->insert(std::make_pair(key, std::move(vec)) );
      }
      
    }else{
      for (size_t cc=0; cc < vv.size(); ++cc){
	auto key = Cycle::hashCode(getKpix(kk.at(cc)),
				   getChannel(kk.at(cc)),
				   bb);
	auto val = vv.at(cc);
	if (target->count(key))
	  target->at(key).push_back(std::move(val));
	else{
	  std::vector<uint16_t> vec;
	  vec.push_back(val);
	  target->insert(std::make_pair(key, std::move(vec)));
	}
      }
      
    }
    
  }
  // just append the vector for each channel
  
}

void rawData::doRmPed(bool rmMAD0){

  Cycle::CalPed();
  size_t                  cySize = m_v_cycles.size();  //
  uint                 currPos=0.0;   //
  uint                   currPct;
  
  // remove ped:
  for (auto &cy : m_v_cycles){
    if (!cy.m_has_adc){
      continue; // skip empty cycles
    }
    
    cy.RemovePed(Cycle::s_ped_adc);  

    currPos++;
    currPct = (uint) (( (double)currPos/(double)cySize )*100.0);
    cout << "\rProcessing cycles to remove pedestal: "
	 << currPct << "%  " << flush;
  }
  cout << "\n";
  
}

void Cycle::RemovePed(std::unordered_map<uint, uint16_t> &ped_adc){
  for (uint b=0; b < m_nbuckets; b++){ // loop over buckets
      auto vv = vadc(b);
      auto kk = hashkeys(b);
      for ( size_t cc=0; cc<vv.size(); ++cc){
    	auto val = vv.at(cc);
	auto key = (getKpix(kk.at(cc)),
		    getChannel(kk.at(cc)),
		    b);
	auto ped = ped_adc.at(key);
      
    	vv.at(cc) = val - ped;
	
      }
    }
  
}

void Cycle::RemoveCM_CalFc(std::unordered_map<uint, double> &slopes, bool remove_adc){

  for (uint b = 0; b< m_nbuckets; b++){

  }
  
}
