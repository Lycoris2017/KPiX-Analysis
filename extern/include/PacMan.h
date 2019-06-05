#ifndef __PACMAN_H__
#define __PACMAN_H__
#include "clustr.h"

using namespace std;


class PacMan
{
	public:
	// Constructor
		PacMan();	
	// Deconstructor
		~PacMan();
		
		double getClusterCoG();
		double getClusterCoGSoN();
		double getClusterCharge();
		double getClusterSignificance();
		int getClusterElementssize();
		
		clustr getCluster();
		void Eater(clustr&, int, int);
		
	private:
		clustr Cluster;

};

#endif
