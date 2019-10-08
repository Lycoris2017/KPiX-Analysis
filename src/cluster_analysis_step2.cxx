//-----------------------------------------------------------------------------
// File          : analysisExternal.cxx
// Author        : Uwe Kraemer (orig. Ryan Herbst) <uwe.kraemer@desy.de>
// Created       : 06/28/2017
// Project       : KPiX Analysis
//-----------------------------------------------------------------------------
// Description :
// Analysis of external triggering KPiX Data.
//-----------------------------------------------------------------------------
// Copyright (c) 2012 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 05/30/2012: created
// 06/28/2017: large scale rewrite of original calibrationFitter.cpp
// 22/03/2018: clean up ecal plots and add strip plots by <mengqing.wu@desy.de>
// 24/09/2018: Branch off of original analysis.cxx to better focus on external trigger data.
//-----------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TTree.h>
#include <TROOT.h>
#include <TPolyLine3D.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TVector.h>
#include <TKey.h>
#include <TClass.h>
#include <stdarg.h>
#include <KpixEvent.h>
#include <KpixSample.h>
#include <Data.h>
#include <DataRead.h>
#include <math.h>
#include <fstream>
#include <XmlVariables.h>
#include <string.h>
#include <algorithm>
#include <stdlib.h>

#include "kpixmap.h"
#include "kpix_left_and_right.h"
#include "testbeam201907_noise_mask.h"
#include "clustr.h"
#include "PacMan.h"
#include "TBFunctions.h"
using namespace std;


//////////////////////////////////////////
// Global Variables
//////////////////////////////////////////
struct tree_cluster_input
{
	double CoG ;
	double Charge;
	double Significance2; //Adding noises in quadrature and later dividing chargesum by this noise sum
	double Sigma;
	int Size;
};

//////////////////////////////////////////
// Functions
//////////////////////////////////////////


//////////////////////////////////////////
// Begin of main analysis
//////////////////////////////////////////
int main ( int argc, char **argv )
{
	TH1::SetDefaultSumw2();

	const unsigned int		n_buckets = 4;
	const unsigned int		n_kpix = 32;//24;
	const unsigned int		n_sensor = n_kpix/2;//24;
//	const unsigned int		n_blocks = 32;
	const unsigned int		n_channels = 1024;
	const unsigned int		n_BCC = 8192;
//	const unsigned int		n_strips = 1840;

	bool					sensorFound[n_sensor] = {false};
	stringstream			FolderName;

	TFile					*rFile;
	int eventNumber, sensor, eventTime;
	tree_cluster_input clusters;
	stringstream           tmp;
	string                 outRoot;

	TH1F 					*cluster_position_y[n_sensor][2];
	//TH1F 					*cluster_position_x[n_sensor][3];
	TH1F 					*num_clusters[n_sensor][2]; // only 2 because it makes no sense to ask for the number of clusters in the seed case cause it will always be 1.
	TH1F 					*cluster_charge[n_sensor][2];
	TH1F 					*cluster_significance2[n_sensor][2];
	TH1F 					*cluster_size[n_sensor][2];
	TH1F 					*cluster_sigma[n_sensor][2];


	string InName = argv[1];

	size_t name_start  = InName.find("/Run") + 1;
	size_t name_length = InName.find(".Tree") - name_start;

	InName = InName.substr(name_start, name_length);

	cout << "Name of output file is " <<  InName << endl;
	tmp.str("");
	tmp << InName << ".cluster_analysed.root";
	outRoot = tmp.str();

	rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
	rFile->cd(); // move into root folder base

	//////////////////////////////////////////
	// Read Data
	//////////////////////////////////////////
	TFile *clusterFile = TFile::Open(argv[1]);
	TTree* clusterTree = (TTree*)clusterFile->Get("clusterTree");

	clusterTree->SetBranchAddress("CoG", &clusters.CoG);
	clusterTree->SetBranchAddress("Size", &clusters.Size);
	clusterTree->SetBranchAddress("Significance2", &clusters.Significance2);
	clusterTree->SetBranchAddress("Sigma", &clusters.Sigma);
	clusterTree->SetBranchAddress("Charge", &clusters.Charge);
	clusterTree->SetBranchAddress("eventNumber", &eventNumber);
	clusterTree->SetBranchAddress("eventTime", &eventTime);
	clusterTree->SetBranchAddress("sensor", &sensor);

	long int NEntries = clusterTree->GetEntries();
	for (long int i = 0; i < NEntries; ++i)
	{
		clusterTree->GetEntry(i);
		cout << "Sensor " << sensor << endl;
		sensorFound[sensor]          = true;

	}

	for (sensor = 0; sensor < n_sensor; sensor++) //looping through all possible kpix
	{
		if (sensorFound[sensor])
		{
			rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
			FolderName.str("");
			FolderName << "Sensor_" << sensor;
			rFile->mkdir(FolderName.str().c_str());
			TDirectory *sensor_folder = rFile->GetDirectory(FolderName.str().c_str());
			sensor_folder->cd();


			tmp.str("");
			tmp << "cluster_position_y_all_sensor" << sensor << "_b0";
			cluster_position_y[sensor][0] = new TH1F(tmp.str().c_str(), "cluster position y1; #mum; #Entries", 1840,-46000, 46000);
			tmp.str("");
			tmp << "cluster_charge_all_sensor" << sensor << "_b0";
			cluster_charge[sensor][0] = new TH1F(tmp.str().c_str(), "cluster charge1; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance2_all_sensor" << sensor << "_b0";
			cluster_significance2[sensor][0] = new TH1F(tmp.str().c_str(), "cluster significance2_1; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_sigma_all_sensor" << sensor << "_b0";
			cluster_sigma[sensor][0] = new TH1F(tmp.str().c_str(), "cluster sigma1; #sigma; #Entries", 200,-0.5, 4.95);
			tmp.str("");
			tmp << "cluster_size_all_sensor" << sensor << "_b0";
			cluster_size[sensor][0] = new TH1F(tmp.str().c_str(), "cluster size1; Size; #Entries", 10,-0.5, 9.5);
			tmp.str("");
			tmp << "clusters_all_sensor" << sensor << "_b0";
			num_clusters[sensor][0] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);

			tmp.str("");
			tmp << "cluster_position_y_CUTS_sensor" << sensor << "_b0";
			cluster_position_y[sensor][1] = new TH1F(tmp.str().c_str(), "cluster position y2; #mum; #Entries", 1840,-46000, 46000);
			tmp.str("");
			tmp << "cluster_charge_CUTS_sensor" << sensor << "_b0";
			cluster_charge[sensor][1] = new TH1F(tmp.str().c_str(), "cluster charge2; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance2_CUTS_sensor" << sensor << "_b0";
			cluster_significance2[sensor][1] = new TH1F(tmp.str().c_str(), "cluster significance2_2; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_sigma_CUTS_sensor" << sensor << "_b0";
			cluster_sigma[sensor][1] = new TH1F(tmp.str().c_str(), "cluster sigma2; #sigma; #Entries", 200,-0.5, 4.95);
			tmp.str("");
			tmp << "cluster_size_CUTS_sensor" << sensor << "_b0";
			cluster_size[sensor][1] = new TH1F(tmp.str().c_str(), "cluster size2; Size; #Entries", 10,-0.5, 9.5);
			tmp.str("");
			tmp << "clusters_CUTS_sensor" << sensor << "_b0";
			num_clusters[sensor][1] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);
		}
	}



	gROOT->GetListOfFiles()->Remove(clusterFile);
	clusterFile->Close();

	rFile->Write();
	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
	rFile->Close();
	
	cout << endl;
	cout << "Writing root plots to " << outRoot << endl;
	cout << endl;

	
	return(0);
}
