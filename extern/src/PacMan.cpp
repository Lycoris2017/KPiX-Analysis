

  // ,##.                   ,==.
 //,#    #.                 \ o ',
//#        #     _     _     \    \
//#        #    (_)   (_)    /    ; 
 //`#    #'                 /   .'  
   //`##'                   "=="


#include "PacMan.h"

using namespace std;

PacMan::PacMan()
{
	Cluster.CoG = 0;
	Cluster.Charge = 0;
}

PacMan::~PacMan() { }




double PacMan::getClusterCoG()
{
	return Cluster.CoG;
}

double PacMan::getClusterCoGSoN()
{
	return Cluster.CoGSoN;
}

double PacMan::getClusterSigma()
{
	return Cluster.Sigma;
}

double PacMan::getClusterCharge()
{
	return Cluster.Charge;
}

double PacMan::getClusterSignificance()
{
	return Cluster.Significance;
}

double PacMan::getClusterSignificance2()
{
	return Cluster.Significance2;
}

int PacMan::getClusterElementssize()
{
	return Cluster.Elements.size(); //divided by two because we have 1 entry for signal and 1 entry for noise
}

clustr PacMan::getCluster()
{
	return Cluster;
}

void PacMan::Eater(clustr& PACMAN, int element, int oldelement, double old_significance) 
{
	if (PACMAN.Elements.count(element) != 0 ) // I always add the starting element and only add the new element if it exists.
	{
		//cout << "DEBUG PACMAN: 1 " << endl;
		double significance = PACMAN.Elements.at(element)/PACMAN.Noise.at(element);
		//cout << "DEBUG PACMAN: 1.1 " << endl;
		if (significance <= old_significance)
		{
			if (PACMAN.Elements.count(element) != 0) //if the current element exists
			{
			
				std::map<int,double>::iterator begin = PACMAN.Elements.begin();
				std::map<int,double>::iterator end = PACMAN.Elements.end();
				//cout << "DEBUG PACMAN: 2 " << endl;
				if (element == begin->first && PACMAN.Elements.size() != 1) // if we are at element 0 of the list but the size is larger than 1
				{
					//cout << "DEBUG PACMAN: 2.1 " << endl;
					if (oldelement - 1 != element) // if we did not come to this point by moving down
					{
						//cout << "Up" << endl;
						PacMan::Eater(PACMAN, element+1, element, significance); //start moving up
					}
				}
				
				if ( element != begin->first && element != end->first) // if we are neither at the first nor at the last element/strip of the list 
				{
					//cout << "DEBUG PACMAN: 2.2 " << endl;
					if (oldelement - 1 == element) // if we came here by moving down
					{
						//cout << "DEBUG PACMAN: 2.2.1 " << endl;
						//cout << "Down" << endl;
						PacMan::Eater(PACMAN, element-1, element, significance); // keep moving down
					}
					if (oldelement + 1 == element) // if we came here by moving up
					{
						//cout << "DEBUG PACMAN: 2.2.2 " << endl;
						//cout << "Up" << endl;
						PacMan::Eater(PACMAN, element+1, element, significance); // keep moving up
					}
					else if (oldelement == 9999)  // starting point for program
					{
						//cout << "DEBUG PACMAN: 2.2.3 " << endl;
						//cout << "Start" << endl;
						PacMan::Eater(PACMAN, element-1, element, significance); // move down
						PacMan::Eater(PACMAN, element+1, element, significance); // move up
					}
				}
				if (element == end->first && PACMAN.Elements.size() != 1) // if we are at the last element of the list but the size is larger than 1
				{
					//cout << "DEBUG PACMAN: 2.3 " << endl;
					//cout << "Down" << endl;
					PacMan::Eater(PACMAN, element-1, element, significance); //start moving down
				}
				
				
				
			}
		
			//cout << "DEBUG PACMAN: 4 " << element << " " << PACMAN.Elements.size() << " " << PACMAN.Elements.count(element) << " " << oldelement << endl;
			Cluster.Elements.insert(std::pair<int,double>(element,PACMAN.Elements.at(element)));
			Cluster.Noise.insert(std::pair<int,double>(element,PACMAN.Noise.at(element)));
			PACMAN.Erase(element);
			//cout << "DEBUG PACMAN: 5 " << endl;
		}
	
		else
		{
			if (oldelement == 9999)
			{
				cout << "charge: " << PACMAN.Elements.at(element) << " | noise: " << PACMAN.Noise.at(element) << " | element: " << element << endl;
				cout << "Error: Starting element has significance higher than 99999" << endl;
			}
			//cout << "Just checking how often it happens that a new strip has higher significance than the previous" << endl;
		}
	}
	if (oldelement == 9999) 
	{
		//cout << "End" << endl;
		Cluster.SetParameters();
	}
	
}
