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
		double getClusterCharge();
		int getElementssize();
		
		clustr getCluster();
		void Eater(clustr&, int, int);
		void Eater_w_Noise(clustr&, int, int);
		
	private:
		clustr Cluster;

};

#endif
