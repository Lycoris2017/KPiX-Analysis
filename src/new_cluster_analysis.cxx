//-----------------------------------------------------------------------------
// Mengqing Wu <mengqing.wu@desy.de>
// @ 2019-12-06
// Target: analysis template using new lycoris database classes 
//-----------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>

#include "rawData.h"
#include "DataRead.h"

#include "kpix_left_and_right.h"
#include "clustr.h"
#include "cluster.h"
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
	
	// gROOT->ProcessLine(".L /home/lycoris-dev/KPiX-Analysis/core/include/cluster.h");
    if (argc != 3){
        printf("[Usage] ./new_analysis [input.dat] [calibration.ymlCalib.root]\n");
		return 0;
	}
	
    printf("[Info] You chose file %s\n", argv[1]);
	rawData db;
    db.setMaxCycles(35000);
    db.setNBuckets(1);
	db.loadFile(argv[1]);
	
	cout<< "[dev] How many cycles? "  << db.getNCycles() << std::endl;
	//db.loadCalib("/home/lycoris-dev/workspace/kpix-analysis/data/calib_HG_20190710T24.csv");

    db.loadCalibTree(argv[2]); // reads data from calibration TTree produced by ymlCalibrationFitter

	db.doRmPedCM();
	//	db.loadGeo("/home/lycoris-dev/workspace/kpix-analysis/data/plane_Geo_default.txt");
    db.loadGeo("/home/lycoris-dev/KPiX-Analysis/data/plane_Geo_default.txt"); // KPiX slot to plane assignment
	uint b = 1;
	Cycle::CalNoise(b);

	auto noisemap = Cycle::getNoise();
	auto timenoisemap = Cycle::getTimeNoise();
    std::string OutRoot = argv[1];
    OutRoot = OutRoot + ".new_cluster.root";

    std::string OutGBL = argv[1];
    OutGBL = OutGBL + ".new_GBL_input.txt";

    TFile *fout = new TFile(OutRoot.c_str(), "recreate");
	fout->cd();
	ofstream claus_file(OutGBL);
	uint header = 1;
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
	std::vector<double> vector_charge, vector_sigma, vector_size, vector_pos, vector_noise, strip_pos, strip_signi, strip_charge;;
	std::vector<int> vector_sensor, vector_ID, vector_kpix, vector_channel, vector_strip, strip_sensor, strip_ID;
    int kpix, channel, eventnumber;
    int sensor, strip, time;
	const uint n_kpix = 24;

    TTree* vectorTree = new TTree("vectorCluster", "vector cluster tree");
    vectorTree->Branch("eventnumber", &eventnumber, "eventnumber/I");
    vectorTree->Branch("sensor", &vector_sensor);
    vectorTree->Branch("charge", &vector_charge);
    vectorTree->Branch("sigma", &vector_sigma);
    vectorTree->Branch("size", &vector_size);
    vectorTree->Branch("time", &time);
    vectorTree->Branch("pos", &vector_pos);
    vectorTree->Branch("ID", &vector_ID);

    TTree* stripTree = new TTree("vectorClusterStrip", "vector strip tree");
    stripTree->Branch("sensor", &strip_sensor);
    stripTree->Branch("charge", &strip_charge);
    stripTree->Branch("signi", &strip_signi);
    stripTree->Branch("pos", &strip_pos);
    stripTree->Branch("ID", &strip_ID);


    TTree* noiseTree = new TTree("noiseTree", "noise tree");
    noiseTree->Branch("eventnumber", &eventnumber, "eventnumber/I");
    noiseTree->Branch("noise", &vector_noise);
    noiseTree->Branch("kpix", &vector_kpix);
    noiseTree->Branch("strip", &vector_strip);
    noiseTree->Branch("channel", &vector_channel);
    noiseTree->Branch("time", &time);


    uint64_t prev_runtime = 0;
    uint ClusterID = 0;
	for (const auto &ev: db.getCycles()){
		if (!ev.m_has_fc) continue;
		eventnumber = ev.eventNumber();
		// fc indexed by bucket+kpix+channel
		std::map<int, double> cluster_Events_after_cut[n_kpix/2];
		std::map<int, double> cluster_Noise_after_cut[n_kpix/2];
		
		//!- 1) fill the input maps
		for (auto &fc: ev.m_m_fc){
			auto key = fc.first; 
			if (Cycle::getBucket(key)!=0) continue;
			
			kpix = Cycle::getKpix(key);
			channel = Cycle::getChannel(key);
			if (kpix%2 ==0 ) strip = kpix2strip_left.at(channel);
			else strip= kpix2strip_right.at(channel);
            time = Cycle::getTime(fc.first);
			sensor = db.getPlane(kpix);

            charge = fc.second;
			noise = noisemap.at(key);
//            cout << "Noise debug: " << noise << endl;
//            cout << "Key debug: " << key << endl;
			// if (eventnumber == 9988 && sensor == 0 && strip ==940){
			// 	printf("debug: ev 9988 sensor 0 strip 940 charge %.4f, noise %.4f\n",
			// 	       charge, noise);
			// }
				
			//!- If no pre-selection on the inputs for cluster, the PacMan algorithm will not erase the input element, thus jump into an infinitive while loop!
            if (charge > 2*noise && strip != 9999 && noise_mask[sensor].at(strip) == 1 ){
				cluster_Events_after_cut[sensor].emplace(strip, charge);
				cluster_Noise_after_cut[sensor].emplace(strip, noise);
			}
            vector_kpix.push_back(kpix);
            vector_channel.push_back(channel);
            vector_strip.push_back(strip);
            vector_noise.push_back(noise);
		}

		// for (sensor =0; sensor< n_kpix/2; sensor++){
		// 	printf("debug: ev %d sensor %d has cluster inputs %d elements, %d noise\n",
		// 	       eventnumber, sensor, cluster_Events_after_cut[sensor].size(),
		// 	       cluster_Noise_after_cut[sensor].size());
		// }

		//!- 2) do clustering sensor by sensor
		for (sensor = 0; sensor < n_kpix/2; sensor++){
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
                NomNom.Eater(Input, Input.MaxSoN(), 9999, 99999, ClusterID);
                ClusterID++;
                clustr Cluster = NomNom.getCluster();
//                cout << "Test" << endl;
                vector_charge.push_back(NomNom.getClusterCharge());
                vector_sigma.push_back(NomNom.getClusterSignificance2());
                vector_size.push_back(NomNom.getClusterElementssize());
                vector_pos.push_back(yParameterSensor(NomNom.getClusterCoG(), sensor));
                vector_ID.push_back(NomNom.getClusterID());

                for (auto const& i : Cluster.Elements)
                {
                    strip_pos.push_back(yParameterSensor(i.first,sensor));
                    strip_charge.push_back(i.second);
                    strip_signi.push_back(i.second/Cluster.Noise.at(i.first));
                    strip_sensor.push_back(sensor);
                    strip_ID.push_back(Cluster.ID);
                }

                stripTree->Fill();
                strip_pos.clear();
                strip_charge.clear();
                strip_signi.clear();
                strip_sensor.clear();
                strip_ID.clear();

                // claus_file:
                if (header == 1){
	                header = 0;
                    claus_file <<"Event Number,Layer,position,Significance,,Size,Charge,runtime,runtime_ns,trigN,ID" << endl;
                }
                //! if you are looking at bucket 0:
                auto trigger0 = ev.m_v_exttrigs.at(0);
                claus_file << setw(5) << eventnumber  << ", "
                           << setw(1) << sensor  << ", "
                           << setw(7) << yParameterSensor(NomNom.getClusterCoG(), sensor) << ", "
                           << setw(7) << NomNom.getClusterSignificance2() << ", "
                           << setw(2) << NomNom.getClusterElementssize() << ", "
                           << setw(7) << NomNom.getClusterCharge() << ", "
                           << setw(13) << std::scientific << std::setprecision(13) << trigger0.runtime << ", "
                           << setw(13) << std::scientific << std::setprecision(13) << (trigger0.runtime)*5.0 << ", "
                           << setw(13) << trigger0.triggerid <<", "
                           << setw(7) << NomNom.getClusterID()
                           << endl;
                
				num_of_clusters++;
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
            vector_sensor.push_back(sensor);
			
		}
        vectorTree->Fill();
        noiseTree->Fill();

        vector_charge.clear();
        vector_sigma.clear();
        vector_size.clear();
        vector_pos.clear();
        vector_sensor.clear();
        vector_ID.clear();

        vector_kpix.clear();
        vector_channel.clear();
        vector_strip.clear();
        vector_noise.clear();

	}


	fout->Write();
	fout->Close();
	printf("[INFO] Clustering ended...\n");
	/*-----------End of Cluster------------*/
    printf("File saved to %s \n", OutRoot.c_str());
    printf("GBL input file saved to %s \n", OutGBL.c_str());
	return 1;
	
}
