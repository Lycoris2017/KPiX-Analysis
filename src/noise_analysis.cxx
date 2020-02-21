//-----------------------------------------------------------------------------
// Mengqing Wu <mengqing.wu@desy.de>
// @ 2019-12-06
// Target: analysis template using new lycoris database classes 
//-----------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string.h>

#include "rawData.h"
#include "DataRead.h"

#include "kpix_left_and_right.h"
#include "clustr.h"
#include "PacMan.h"
#include "TBFunctions.h"
#include "testbeam201907_noise_mask.h"

#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"

using namespace std;
using namespace Lycoris;

// Process the data
int main ( int argc, char **argv ) {
	auto kpix2strip_left = kpix_left();
	auto kpix2strip_right = kpix_right();
	
    gROOT->ProcessLine(".L /home/lycoris-dev/KPiX-Analysis/core/include/cluster.h");
	if (argc == 1){
        printf("[Usage] ./new_cluster_analysis [input.dat] [calib.root] \n");
		return 0;
	}
	
	printf("[Info] You choose file %s\n", argv[1]);
	rawData db;
	//db.setMaxCycles(1000);
	db.setNBuckets(1);
	db.loadFile(argv[1]);
	
	cout<< "[dev] How many cycles? "  << db.getNCycles() << std::endl;
	//db.loadCalib("/home/lycoris-dev/workspace/kpix-analysis/data/calib_HG_20190710T24.csv");
    if (argc != 3){
        db.loadCalib("/home/lycoris-dev/workspace/kpix-analysis/data/HG_slopes_D.root");
    }
    else{
        db.loadCalibTree(argv[2]);
    }


	db.doRmPedCM();
	db.loadGeo("/home/lycoris-dev/workspace/kpix-analysis/data/plane_Geo_default.txt");
    uint b = 4;
	Cycle::CalNoise(b);
    auto noiseTimemap = Cycle::getTimeNoise();
    auto noisemap = Cycle::getNoise();
    std::string OutRoot = argv[1];
    OutRoot = OutRoot + ".noise.root";
//    std::string OutGBL = argv[1];
//    OutGBL = OutGBL + ".new_GBL_input.txt";

    TFile *fout = new TFile(OutRoot.c_str(), "recreate");
	fout->cd();
	// TH2F *cluster_correlation_s0_s1 = new TH2F("cluster_correlation_s0_s1",
	//                                            "strip correlation; sensor0 [#mum]; sensor1[#mum]",
	//                                            230, -46000.5, 45999.5, 230,-46000, 46000);

	/*-----------Start of Cluster------------*/
    printf("[INFO] Noise writedown started...\n");
	
	
    std::vector<double> noise, noise2;
    int eventnumber;
    std::vector<int> kpix, channel, strip, bucket, time;
	const uint n_kpix = 24;

    TTree* noiseTimeTree = new TTree("noise_time_tree", "noise_time_tree");
    noiseTimeTree->Branch("eventnumber", &eventnumber, "eventnumber/I");
    noiseTimeTree->Branch("kpix", &kpix);
    noiseTimeTree->Branch("channel", &channel);
    noiseTimeTree->Branch("strip", &strip);
    noiseTimeTree->Branch("bucket", &bucket);
    noiseTimeTree->Branch("time", &time);
    noiseTimeTree->Branch("noise", &noise);

    TTree* noiseTree = new TTree("noise_tree", "noise_tree");
    noiseTree->Branch("eventnumber", &eventnumber, "eventnumber/I");
    noiseTree->Branch("kpix", &kpix);
    noiseTree->Branch("channel", &channel);
    noiseTree->Branch("strip", &strip);
    noiseTree->Branch("bucket", &bucket);
    noiseTree->Branch("noise", &noise2);

	for (const auto &ev: db.getCycles()){
		if (!ev.m_has_fc) continue;
		eventnumber = ev.eventNumber();
		// fc indexed by bucket+kpix+channel
		for (auto &fc: ev.m_m_fc){
			auto key = fc.first;
            auto k =Cycle::getKpix(key);
            auto c = Cycle::getChannel(key);
            auto rTkey = Cycle::rmTime(fc.first);

            if (k%2 ==0 ) strip.push_back(kpix2strip_left.at(c));
            else strip.push_back(kpix2strip_right.at(c));
            kpix.push_back(Cycle::getKpix(key));
            channel.push_back(Cycle::getChannel(key));
            bucket.push_back(Cycle::getBucket(key));
            time.push_back(Cycle::getTime(key));
            noise.push_back(noiseTimemap.at(key));
            noise2.push_back(noisemap.at(rTkey));


		}
        noiseTree->Fill();
        noiseTimeTree->Fill();

        kpix.clear();
        channel.clear();
        bucket.clear();
        time.clear();
        noise.clear();
        noise2.clear();
        strip.clear();
	}

	fout->Write();
	fout->Close();
    printf("[INFO] Noise writedown ended...\n");
	/*-----------End of Cluster------------*/
    printf("File saved to %s \n", OutRoot.c_str());
//    printf("GBL input file saved to %s \n", OutGBL.c_str());
	return 1;
	
}
