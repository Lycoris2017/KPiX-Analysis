#ifndef __PACMAN_H__
#define __PACMAN_H__
#include "clustr.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_RESET   "\x1b[0m"

using namespace std;


class PacMan
{
	public:
	// Constructor
		PacMan();	
	// Deconstructor
		~PacMan();
        uint getClusterID();
		double getClusterCoG();
		double getClusterSigma();
		double getClusterCoGSoN();
		double getClusterCharge();
		double getClusterSignificance();
		double getClusterSignificance2();
		int getClusterElementssize();
		
		clustr getCluster();
        void Eater(clustr&, int, int, double, uint);
        void Eater(clustr&, int, int, double);
		
	private:
		clustr Cluster;

};

#endif
