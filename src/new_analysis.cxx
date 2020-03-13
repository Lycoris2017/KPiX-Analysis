//-----------------------------------------------------------------------------
// Mengqing Wu <mengqing.wu@desy.de>
// @ 2019-10-16
// Target: analysis template using new lycoris database classes
//-----------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string.h>

#include "rawData.h"
#include "DataRead.h"

#include "kpix_left_and_right.h"

#include "TTree.h"
#include "TFile.h"

using namespace std;
using namespace Lycoris;

// Process the data
int main ( int argc, char **argv ) {
	auto kpix2strip_left = kpix_left();
	auto kpix2strip_right = kpix_right();
	
	
    if (argc != 3){
        printf("[Usage] ./new_analysis [input.dat] [calibration.ymlCalib.root]\n");
		return 0;
	}
	
	printf("[Info] You choose file %s\n", argv[1]);
	rawData db;
	//db.setMaxCycles(100);
	db.setNBuckets(1);
	db.loadFile(argv[1]);
	
	cout<< "[dev] How many cycles? "  << db.getNCycles() << std::endl;
	//db.loadCalib("/home/lycoris-dev/workspace/kpix-analysis/data/calib_HG_20190710T24.csv");

    db.loadCalibTree(argv[2]);


	//return(0);
	db.doRmPedCM();
	
	uint b = 1;
	Cycle::CalNoise(b);
	auto noisemap = Cycle::getNoise();
	auto pedestal = Cycle::s_ped_med_adc;
	auto mads     = Cycle::s_ped_mad_adc;
    auto calibs = db.getCalibs();
	if (noisemap.empty()) return 0;
	
	
	TFile *fout = new TFile("test.root", "recreate");
	fout->cd();
	double noise, ped, mad;
	int kpix, bucket, strip, channel;
	TTree* test = new TTree("ttree", "test tree");
	test->Branch("kpix",   &kpix,   "kpix/I");
	test->Branch("strip",  &strip,  "strip/I");
	test->Branch("bucket", &bucket, "bucket/I");
	test->Branch("channel", &channel, "channel/I");
	test->Branch("ped",    &ped,    "ped/D"); //fc
	test->Branch("mad",    &mad,    "mad/D"); //fc
	test->Branch("noise",  &noise,  "noise/D"); //fc
	// loop to fill the tree:
	for (const auto &a: pedestal){
		
		auto key = a.first;
		bucket = Cycle::getBucket(key);
		if (bucket!=0) continue;
		ped = (double)pedestal.at(key);
		mad = (double)mads.at(key);
		if (noisemap.count(key)) noise = noisemap.at(key);
		else noise = 999;
		
		kpix = Cycle::getKpix(key);
		channel = Cycle::getChannel(key);
		if (kpix%2) strip = kpix2strip_left.at(channel);
		else strip = kpix2strip_right.at(channel);
        if (calibs.at(key).first==0) {
			std::cout <<"k " << kpix << " c "<< channel
                      <<" s " << calibs.at(key).first << std::endl;
			//continue;
		}
		
		//ped = ped/slopes.at(key);
        mad = mad/calibs.at(key).first;
		
		test->Fill();
	}

	double cm_med, charge_fc;
	int kpix1, bucket1, channel1, strip1, eventnumber;
	TTree* cycles = new TTree("cycles","cycle tree");
	cycles->Branch("charge_fc", &charge_fc, "charge_fc/D");
	cycles->Branch("kpix",   &kpix1,   "kpix1/I");
	cycles->Branch("strip",  &strip1,  "strip1/I");
	cycles->Branch("channel",  &channel1,  "channel1/I");
	cycles->Branch("bucket", &bucket1, "bucket1/I");
	cycles->Branch("eventnumber", &eventnumber, "eventnumber/I");
	cycles->Branch("cm_med", &cm_med, "cm_med/D");

	for (const auto &ev: db.getCycles()){
		if (!ev.m_has_fc) continue;
		bucket1=0;
		eventnumber = ev.m_cyclenumber;
		for (auto &fc: ev.m_m_fc){
			auto key = fc.first;
			kpix1 = Cycle::getKpix(key);
			channel1 = Cycle::getChannel(key);
			if (kpix%2) strip = kpix2strip_left.at(channel1);
			else strip1 = kpix2strip_right.at(channel1);

			charge_fc = fc.second;
			cm_med = ev.m_m_cm_noise.at(kpix1);
			
			cycles->Fill();
		}
	}

	
	fout->Write();
	fout->Close();
	return 1;
	
}
