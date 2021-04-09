#include <KpixEvent.h>
#include <KpixSample.h>
#include <fstream>
#include <cstring>
#include <chrono>
#include <sstream>
#include <cassert> // assert

#include "TFile.h"
#include "TKey.h"
#include "TList.h"
#include "TH1.h"
#include "TObject.h"
#include "TTree.h"

#include "rawData.h"
#include "TBFunctions.h"
#include "Data.h"
#include "DataRead.h"


using namespace std;
using namespace Lycoris;

/* !ALL static variables need declaration before using! */
//std::unordered_map< uint, std::vector<double> > Cycle::s_buf_fc;
//std::unordered_map< uint, std::vector<double>> Cycle::s_buf_time_fc;
std::unordered_map< uint, std::unordered_map<uint, std::vector<double>> > Cycle::s_buf_time_fc;
std::unordered_map< uint, double > Cycle::s_noise_fc;
//std::unordered_map< uint, double > Cycle::s_noise_time_fc;
std::unordered_map< uint, std::unordered_map<uint, double> > Cycle::s_noise_time_fc;

std::unordered_map< uint, std::vector<int> > Cycle::s_buf_adc;
std::unordered_map< uint, double > Cycle::s_ped_med_adc;
std::unordered_map< uint, double > Cycle::s_ped_mad_adc;

Cycle::Cycle(KpixEvent &event,
             uint &ntrig_ext,
             uint nbuckets,
             uint begin_ch,
             uint end_ch,
             bool isold ){

	// Declare member variables first
	m_nbuckets = nbuckets;
	m_has_adc = false;
	m_has_fc  = false;
	m_cyclenumber = event.eventNumber();
    if (isold) m_ts64 = event.timestamp();
    else       m_ts64 = event.runtime();
	
	
	KpixSample *sample;   //
	for (uint ev=0; ev< event.count(); ev++){
		
		sample = event.sample(ev);
		
		uint kpix, channel, bucket, value;
		uint bunchClk, subCount, tstamp;
		kpix    = sample->getKpixAddress();
		channel = sample->getKpixChannel();
		bucket  = sample->getKpixBucket();
		value   = sample->getSampleValue();
		tstamp  = sample->getSampleTime(); // local time in BCC
		
		if (sample->getSampleType() == KpixSample::Data){
			if ( bucket >= m_nbuckets ) continue;
			if ( channel < begin_ch || channel > end_ch)  continue;
			
			uint key = hashCode(kpix, channel, bucket);
			
			m_v_hashkeys.push_back(key);
			m_v_adc.push_back(value);
			m_v_ts.push_back(tstamp);
		}
		
		//! External Trigger Data
		if (sample->getSampleType() == KpixSample::Timestamp){
			double time;
			trigger_t trigger;
			
			//! global trigger counter increment
			ntrig_ext++; 
			trigger.triggerid = ntrig_ext;
			
            if (isold)
                time = tstamp + double(value * 0.125);
            else{
                //tstamp   = sample->getSampleRuntime64(m_ts64);
                subCount = sample->getSubCount();
                bunchClk = sample->getBunchCount();
                time = bunchClk + double(0.125 * subCount);
                trigger.runtime = sample->getSampleRuntime64(m_ts64);
            }
			// TBD: Do nothing now, add it to data structure later.
			trigger.tstamp  = time;
			m_v_exttrigs.push_back(trigger);

			//std::cout<<"[dev] Trigger ts = " << time 
			//	 <<", runtime[ns] = " <<  tstamp*5 << std::endl;
		}
		

    }// finish loop over all data of one cycle

    //printf("adc vec size : %d\n", m_v_adc_b[i].size() );
    m_has_adc = m_has_adc || m_v_adc.size();
    
    // **Exception Handler - Start
    try{
	    if (m_v_hashkeys.size() != m_v_adc.size())
		    throw m_cyclenumber;
    }
    catch (uint e){
	    cout << "ERROR: Cycle Event " << e << " has channel_keys != its value. \n";
    }
    // **Exception Handler - End
}


void rawData::loadCalibTree(const std::string& fname){
    string root = ".root";
    printf(" loadCalibTree ");

    try{
        if (fname.find(root) !=std::string::npos)
            loadRootTree(fname);
        else
            throw std::runtime_error("loadCalibTree: invalid input calib file!\n");
    }
    catch(std::exception &e){
        cout << "Caught exception: " << e.what() << "\n";
    }
}

void rawData::loadRootTree(const std::string & fname){
    printf(" from ROOT... \n");
    TFile *calib = TFile::Open(fname.c_str());
    if (!calib) {
        throw std::runtime_error("loadRootTree: invalid calib root file!\n");
        return;
    }
    printf("Open Calib file : %s\n", fname.c_str());
    m_m_calibs.clear();

    TTree *calib_tree = (TTree*)calib->Get("calibration_tree");
    uint kpix_num, channel_num, bucket_num, range_num;
    double calib_slope, calib_error, calib_pearsson;
    calib_tree->SetBranchAddress("kpix", &kpix_num);
    calib_tree->SetBranchAddress("channel", &channel_num);
    calib_tree->SetBranchAddress("bucket", &bucket_num);
    calib_tree->SetBranchAddress("range", &range_num);
    calib_tree->SetBranchAddress("calib_slope", &calib_slope);
    calib_tree->SetBranchAddress("calib_error", &calib_error);
    calib_tree->SetBranchAddress("calib_pearsson", &calib_pearsson);
    long int nEnTrees = calib_tree->GetEntries();
    for (long int i = 0; i < nEnTrees; ++i){
        calib_tree->GetEntry(i);
        auto index = Cycle::hashCode(kpix_num, channel_num, bucket_num);
        std::pair<double, double> calib_all = std::make_pair(calib_slope, calib_pearsson);
        m_m_calibs.emplace(index, calib_all);
    }
    calib->Close();
}

void rawData::loadCSV(const std::string& fname){
  /* 
     Read Calib from a csv file
     cols needs to be:
     kpix >> channel >> bucket >> slope
  */
	printf(" from CSV...\n");
	m_m_slopes.clear();
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
			//if (bucket==0)
			m_m_slopes.emplace(Cycle::hashCode(kpix,channel,bucket), slope);
		}
	}
  
	cout << " test: "<< m_m_slopes.at(Cycle::hashCode(1, 1000, 0)) << endl;
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

void rawData::loadFile(const std::string& fname){
  DataRead               dataRead;
  off_t                  fileSize;  //
  off_t                  filePos;   //
  uint                   lastPct;
  uint                   currPct;
  KpixEvent              event;    //
  uint                   ncys=0;
  string                 calState;
  uint                   ntrig_ext=0;
  
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
    ncys++;
    if (m_nmax!=0 && ncys > m_nmax) break;

    Cycle cy(event, ntrig_ext, m_nbuckets );
    calState   = dataRead.getYmlStatus("CalState");
    if (calState == "Baseline" || calState == "Idle"){ //Only adding baseline cycles
        Cycle::AddAdcBuf(cy);
        // For pedestal calculation:

        // !!Once after std::move, you should never use this object
        m_v_cycles.push_back(std::move(cy));
    }
    //else cout << "Found " << calState << " cycle." << endl;


    
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

  dataRead.close();
}
//calState   = dataRead.getYmlStatus("CalState");
/* Static */
void Cycle::CalPed(uint nbuckets, bool save_mad){

  if (s_buf_adc.empty()) return ;

  s_ped_med_adc.clear();
  s_ped_mad_adc.clear();
  //  s_v_mad0_chn.clear();

  auto t_start = std::chrono::high_resolution_clock::now();

  for( auto &buf : s_buf_adc){
      uint bucket = getBucket(buf.first);
	  // skip not requested buckets
      if (bucket < nbuckets ) {
		  // buf.first, buf.second -> median
		  double ped = median(&buf.second);
		  double mad = MAD(&buf.second);
		  
		  //printf("debug: buf.second size = %d, median = %d, MAD = %d\n", buf.second.size(), ped, mad);
		  s_ped_med_adc.insert(std::make_pair(buf.first, ped ));
		  if (save_mad) s_ped_mad_adc.insert(std::make_pair(buf.first, mad));
		  //if (mad==0)  s_v_mad0_chn.push_back(buf.first);
	  }
  }
  printf("Pedestal & MAD calculation finished.\t");

  auto t_end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds> (t_end - t_start).count();
  cout << " Time[s]: " << duration
       << endl;

  // after pedestal calculation, clean the buffer:
  Cycle::ResetAdcBuf();

}

/* Static */
void Cycle::AddAdcBuf(Cycle& cy){
  if (!cy.m_has_adc) return;
  //  for (uint bb =0; bb< cy.m_nbuckets; bb++) {
  auto kk = cy.hashkeys();
  auto vv = cy.vadc(); 
  for (size_t cc=0; cc<vv.size(); ++cc){
	  uint key = kk.at(cc); // uint
	  int  val = (int)vv.at(cc); // uint16 
	  Cycle::AddBufferT(Cycle::s_buf_adc, key, val);

  }
}


void rawData::doRmPedCM(bool rmAdc){

	Cycle::CalPed(m_nbuckets, true);
	size_t     cySize = m_v_cycles.size();  //
	uint       currPos=0.0;   //
	uint       currPct;
	
	auto t_start = std::chrono::high_resolution_clock::now();
	
	// remove ped:
	for (auto &cy : m_v_cycles){
		if (!cy.m_has_adc){
			continue; // skip empty cycles
		}
		
		//    cy.RemovePed_CalCM_fC(Cycle::s_ped_med_adc,
		//                          m_m_slopes,
		//                          Cycle::s_ped_mad_adc,
		//                          true);
		cy.RemovePed_CalCM_fC(Cycle::s_ped_med_adc,
		                      m_m_calibs,
		                      Cycle::s_ped_mad_adc,
		                      true);
		cy.RemoveCM();
		
		Cycle::AddFcBuf(cy);
		
		// debug:
		//---> Fancy Bar
		currPos++;
		currPct = (uint) (( (double)currPos/(double)cySize )*100.0);
		cout << "\rProcessing cycles to remove pedestal, adc->fC, ignore MAD==0 : "
		     << currPct << "%  " << flush;
		
		//---> Fancy Bar <---// 
		
  }
	auto t_end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds> (t_end - t_start).count();
	cout << "Finished common mode calculation and pedestal subtraction" << endl;
	cout << " Time[s]: " << duration
	     << endl;
	
	cout << "\n";
	
}

void Cycle::RemovePed_CalCM_fC(std::unordered_map<uint, double> &ped_adc,
                               std::unordered_map<uint, std::pair<double, double>> &calibs,
                               std::unordered_map<uint, double> &ped_mad,
                               bool remove_adc,
                               bool cut_mad0,
                               bool cut_calibs)
{
  /* Only work at bucket 0 */
    //  uint bucket=0;
  // Fill in a buffer to calculate CM, indexed by kpix num:
  std::unordered_map<uint, vector<double>> cm_noise_buf;
  auto target = &m_m_fc;
  auto vv = vadc();
  auto kk = hashkeys();
  for (size_t cc =0; cc<vv.size(); ++cc){
	  
	  uint key  = kk.at(cc);
	  uint adc = vv.at(cc);
	  uint ped = ped_adc.at(key);
	  std::pair<double, double> calib = calibs.at(key);

	  uint kpix = getKpix(key);
	  uint channel = getChannel(key);
	  
	  // ignore mad0 channels
	  if (ped_mad.at(key)==0) continue;
	  // ignore channels with slope ==0 and with a pearson correlation coefficient < 0.8
	  if (calib.first == 0 || calib.second < 0.9) continue;
	  
	  double fc = (double)(1.0*adc - ped) / calib.first;
	  // // debug:
	  // if (m_cyclenumber ==99 && kpix==0){
	  //     printf("cycle 99, kpix %d channel %d with fc %.2f (adc %d, ped %d, slope %.2f)\n",
	  //            kpix, channel, fc, adc, ped, slope);
	  // }
	  
	  target->insert(std::make_pair(key, fc));
	  //    cout << "DEBUG bla" << endl;
	  // if (kpix==1 && getChannel(key) == 73){
	  //     cout << "debug: ev " << m_cyclenumber << " k1 c73: "
	  //          << adc << " adc, "
	  //          << ped << " adc, "
	  //          << slope << " slope, "
	  //          << fc << " fC "
	  //          << endl;
	  // }
	  
	  if (cm_noise_buf.count(kpix))
		  cm_noise_buf.at(kpix).push_back(fc);
	  else{
		  std::vector<double > vec;
		  vec.push_back(fc);
		  cm_noise_buf.insert(std::make_pair(kpix, std::move(vec)));
	  }
//    cout << "DEBUG end bla" << endl;
//    if (m_cyclenumber == 21 && kpix ==1){
//	    cout << "cm input debug: ev 21 kpix 1 channel " << channel << ", "
//	         << adc << " adc, "
//	         << ped << " adc, "
//	         << slope << " slope, "
//	         << fc << " fC "
//	         << endl;
//    }

  }

  
  m_has_fc = !(target->empty());
  if (remove_adc) ResetAdc();

  /* Calculate CM */
  m_m_cm_noise.clear();
  //printf("debug: how many kpix in cm_noise calculation? %d \n", cm_noise_buf.size());
//    cout << "DEBUG end bla bla" << endl;
  //std::map<uint, double> debug;
  for(auto &a:cm_noise_buf ){
      double cm_noise = median(&a.second);

      m_m_cm_noise.insert(std::make_pair(a.first, cm_noise));
      //  debug.insert(std::make_pair(a.first, cm_noise));
 }
  // for (auto &d: debug)
  // 	  cout << "Common modes median of EventNumber " << m_cyclenumber
  // 	       << " kpix "  << d.first
  // 	       << " entry " << d.second << endl; // debug

//  cout << "DEBUG end bla bla blaaaaaa" << endl;
    
}


void Cycle::RemoveCM(){

  if (m_m_fc.empty()) return;

  for ( auto &fc: m_m_fc){
	  auto key = rmBucket(fc.first);
	  uint kpix = getKpix(key);
	  double cm_noise = m_m_cm_noise.at(kpix);
	  fc.second = fc.second - cm_noise;
  }
}

/* Static */
void Cycle::AddFcBuf(Cycle& cy){
	/* Note: every cycle, for each k/ch/b, 
	   charge and BCC tstamp is one-to-one correspondence*/
  if (!cy.m_has_fc) return;
  // for (auto &fc: cy.m_m_fc){
  // 	  //uint key = rmTime(fc.first);
  // 	  // Cycle::AddBufferT( Cycle::s_buf_fc,
  //     //                    key, fc.second); // filling only timestampless key
  //     Cycle::AddBufferT( Cycle::s_buf_time_fc,
  //                        fc.first, fc.second);
  // }

  //! MW: <uint, <uint, vector<T>>>
  auto hwindex = cy.hashkeys();
  auto tstamps = cy.vtstamp();
  for (size_t cc=0; cc<hwindex.size(); ++cc){
	  uint outer = hwindex.at(cc);
	  uint inner = tstamps.at(cc);
	  if (!cy.m_m_fc.count(outer)) continue;
	  double val = cy.m_m_fc.at(outer);
	    
	  if (!s_buf_time_fc.count(outer)){
		  std::unordered_map<uint, vector<double>> innermap;
		  s_buf_time_fc.insert(std::make_pair(outer, std::move(innermap)));
	  }
	    
	  Cycle::AddBufferT( s_buf_time_fc.at(outer),
	                     inner, val);

  }
 
}

/* Static */
void Cycle::CalNoise(uint& nbuckets){
	
  if (s_buf_time_fc.empty()) return;
  
  //  s_noise_fc.clear();
  Cycle::ResetNoise();
  auto t_start = std::chrono::high_resolution_clock::now();

  //! Calculate timeless noise:
//   for (auto &buf: s_buf_fc){
//       uint bucket = getBucket(buf.first);
//       if (bucket >= nbuckets) continue;
      
//       double noise = 1.4826*MAD(&buf.second);
//       s_noise_fc.insert(std::make_pair(buf.first, noise));
      
// //      if (getKpix(buf.first)==1 && getChannel(buf.first)==73){
// //	      cout << "debug: noise cal kpix 1, channel 73 noise " <<  noise << endl;
// //          cout << "Key: " << buf.first << endl;
// //	      for (const auto &input: buf.second)
// //		      cout << " " << input << endl;
// //      }
      
    
//   }
  
  //! Calculate timed noise:
  for (auto &buf: s_buf_time_fc){
	  
      uint bucket = getBucket(buf.first);
      if (bucket >= nbuckets) continue;

      std::vector<double> buf_timeless_fc;
      std::unordered_map<uint, double> noise_time_fc;
      for (auto &buf_time: buf.second){
	      uint ts = buf_time.first;
	      double noise_time = 1.4826*MAD(&buf_time.second);
	      noise_time_fc.insert(std::make_pair(ts, noise_time));

	      //!concatenate all timed fC to one vector
	      buf_timeless_fc.insert(buf_timeless_fc.end(),
	                             buf_time.second.begin(),
	                             buf_time.second.end());
      }
      s_noise_time_fc.insert(std::make_pair(buf.first, std::move(noise_time_fc)));

      //! debug:
      if (buf_timeless_fc.size()==0){
	      std::cout << "[debug] 0 sized timeless fC buffer - kpix: "
	                << getKpix(buf.first)
	                << ", channel " << getChannel(buf.first)
	                << std::endl;
	      continue;
      }
      double noise = 1.4826*MAD(&buf_timeless_fc);
      s_noise_fc.insert(std::make_pair(buf.first, noise));
  }
  
  printf("Noise calculation finished.\t");
  auto t_end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds> (t_end - t_start).count();
  cout << " Time[s]: " << duration
       << endl;

  Cycle::ResetFcBuf();
}
