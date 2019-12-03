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
	
	
	if (argc == 1){
		printf("[Usage] ./new_analysis [input.dat]\n");
		return 0;
	}
	
	printf("[Info] You choose file %s\n", argv[1]);
	rawData db;
	db.setMaxCycles(1000);
	db.setNBuckets(1);
	db.loadFile(argv[1]);
	
	cout<< "[dev] How many cycles? "  << db.getNCycles() << std::endl;
	db.loadCalib("/home/lycoris-dev/workspace/kpix-analysis/data/calib_HG_20190710T24.csv");
	
	db.doRmPedCM();
	
	uint b = 1;
	Cycle::CalNoise(b);
	auto noisemap = Cycle::getNoise();
	auto pedestal = Cycle::s_ped_adc;
	auto mads     = Cycle::s_ped_mad;
	auto slopes = db.getSlopes();
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
		
		key  = Cycle::rmBucket(key);
		kpix = Cycle::getKpix(key);
		channel = Cycle::getChannel(key);
		if (kpix%2) strip = kpix2strip_left.at(channel);
		else strip = kpix2strip_right.at(channel);
		if (slopes.at(key)==0) continue;
		ped = ped/slopes.at(key);
		mad = mad/slopes.at(key);
		
		test->Fill();
	}

	double cm_med;
	int kpix1, bucket1;
	TTree* cycles = new TTree("cycles","cycle tree");
	
	
	
	fout->Write();
	fout->Close();
	return 1;
	
}
