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
	//db.setMaxCycles(1000);
	db.setNBuckets(1);
	db.loadFile(argv[1]);
	
	cout<< "[dev] How many cycles? "  << db.getNCycles() << std::endl;
	//db.loadCalib("/home/lycoris-dev/workspace/kpix-analysis/data/calib_HG_20190710T24.csv");
	db.loadCalib("/home/lycoris-dev/workspace/kpix-analysis/data/HG_slopes_D.root");

	db.doRmPedCM();
	db.loadGeo("/home/lycoris-dev/workspace/kpix-analysis/data/plane_Geo_default.txt");
	uint b = 1;
	Cycle::CalNoise(b);

	auto noisemap = Cycle::getNoise();

	TFile *fout = new TFile("cluster.root", "recreate");
	fout->cd();
	// TH2F *cluster_correlation_s0_s1 = new TH2F("cluster_correlation_s0_s1",
	//                                            "strip correlation; sensor0 [#mum]; sensor1[#mum]",
	//                                            230, -46000.5, 45999.5, 230,-46000, 46000);

	/*-----------Start of Cluster------------*/
	printf("[INFO] Clustering started...\n");
	unordered_map<uint, uint> noise_mask[6];
	noise_mask[0] = noise_sensor_0();
	noise_mask[1] = noise_sensor_1();

	noise_mask[2] = noise_sensor_2();
	noise_mask[3] = noise_sensor_3();
	noise_mask[4] = noise_sensor_4();
	noise_mask[5] = noise_sensor_5();
	
	
	double charge, noise;
	double cluster_charge, cluster_pos, cluster_sigma, cluster_size;
	int kpix, channel, eventnumber;
	int sensor, strip;
	const uint n_kpix = 24;

	TTree* ctree = new TTree("cluster","cluster tree");
	ctree->Branch("eventnumber", &eventnumber, "eventnumber/I");
	ctree->Branch("sensor",      &sensor,      "sensor/I");
	ctree->Branch("cluster_charge", &cluster_charge, "cluster_charge/D");
	ctree->Branch("cluster_pos",    &cluster_pos,    "cluster_pos/D");
	ctree->Branch("cluster_sigma",  &cluster_sigma,  "cluster_sigma/D");
	ctree->Branch("cluster_size",   &cluster_size,   "cluster_size/D");

	for (const auto &ev: db.getCycles()){
		if (!ev.m_has_fc) continue;
		eventnumber = ev.eventNumber();
		// fc indexed by bucket+kpix+channel
		std::map<int, double> cluster_Events_after_cut[n_kpix/2];
		std::map<int, double> cluster_Noise_after_cut[n_kpix/2];
		clustr Cluster[n_kpix/2];
		std::vector<clustr> multi_cluster[8][n_kpix/2];
		
		//!- 1) fill the input maps
		for (auto &fc: ev.m_m_fc){
			auto key = fc.first;
			if (Cycle::getBucket(key)!=0) continue;
			
			kpix = Cycle::getKpix(key);
			channel = Cycle::getChannel(key);

			if (kpix%2 ==0 ) strip = kpix2strip_left.at(channel);
			else strip= kpix2strip_right.at(channel);

			sensor = db.getPlane(kpix);

			charge = fc.second;
			noise = noisemap.at(key);

			// if (eventnumber == 9988 && sensor == 0 && strip ==940){
			// 	printf("debug: ev 9988 sensor 0 strip 940 charge %.4f, noise %.4f\n",
			// 	       charge, noise);
			// }
				
			//!- If no pre-selection on the inputs for cluster, the PacMan algorithm will not erase the input element, thus jump into an infinitive while loop!
			if (charge > 3*noise && strip != 9999 && noise_mask[sensor].at(strip) == 1 ){
				cluster_Events_after_cut[sensor].emplace(strip, charge);
				cluster_Noise_after_cut[sensor].emplace(strip, noise);
			}
		}

		// for (sensor =0; sensor< n_kpix/2; sensor++){
		// 	printf("debug: ev %d sensor %d has cluster inputs %d elements, %d noise\n",
		// 	       eventnumber, sensor, cluster_Events_after_cut[sensor].size(),
		// 	       cluster_Noise_after_cut[sensor].size());
		// }

		//!- 2) do clustering sensor by sensor
		for (sensor =0; sensor< n_kpix/2; sensor++){
			if (cluster_Events_after_cut[sensor].size()==0) continue;
			clustr Input;
			Input.Elements = cluster_Events_after_cut[sensor];
			Input.Noise   = cluster_Noise_after_cut[sensor];
			int num_of_clusters=0;
			//printf("debug- start while\n");
			//- Start clustering
			while(Input.Elements.size()!=0){
				PacMan NomNom;
				//double SoN_order = 0;

				//MaximumSoN[sensor][Input.MaxSoN()]+=weight;
				NomNom.Eater(Input, Input.MaxSoN(), 9999, 99999);

				cluster_charge = NomNom.getClusterCharge();
				cluster_sigma  = NomNom.getClusterSignificance2();
				cluster_size   = NomNom.getClusterElementssize();
				cluster_pos    = yParameterSensor(NomNom.getClusterCoG(), sensor);
				
				// debug to print as claus_file:
				// cout << setw(5) << eventnumber  << ","
				//      << setw(1) << sensor  << ","
				//      << setw(7) << cluster_pos << ","
				//      << setw(7) << cluster_sigma << ","
				//      << setw(2) << cluster_size << ","
				//      << setw(7) << cluster_charge << ","
				//      << endl;
				num_of_clusters++;
				ctree->Fill();
			}
			// printf("debug: ev %d sensor %d has %d elements, %d noise, %d clusters.\n",
			//        eventnumber,
			//        sensor,
			//        cluster_Events_after_cut[sensor].size(),
			//        cluster_Noise_after_cut[sensor].size(),
			//        num_of_clusters);
			
			if (eventnumber==9988 && sensor==0){
				printf("[debug] ev 9988 sensor 0 has input:\n");
				for (const auto &a : cluster_Events_after_cut[sensor])
					printf("\tstrip %d has %.4f charge, %.4f noise\n",
					       a.first,
					       a.second,
					       cluster_Noise_after_cut[sensor].at(a.first));
				
			}
			
		}
	}

	
	fout->Write();
	fout->Close();
	printf("[INFO] Clustering ended...\n");
	/*-----------End of Cluster------------*/
	
	return 1;
	
}
