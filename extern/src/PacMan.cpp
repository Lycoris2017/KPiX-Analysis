

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

double PacMan::getClusterCharge()
{
	return Cluster.Charge;
}

int PacMan::getElementssize()
{
	return Cluster.Elements.size();
}

clustr PacMan::getCluster()
{
	return Cluster;
}

void PacMan::Eater(clustr& PACMAN, int element, int oldelement) 
{
	//cout << "DEBUG PACMAN: 1 " << endl;
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
				PacMan::Eater(PACMAN, element+1, element); //start moving up
			}
		}
		
		if ( element != begin->first && element != end->first) // if we are neither at the first nor at the last element/strip of the list 
		{
			//cout << "DEBUG PACMAN: 2.2 " << endl;
			if (oldelement - 1 == element) // if we came here by moving down
			{
				//cout << "DEBUG PACMAN: 2.2.1 " << endl;
				//cout << "Down" << endl;
				PacMan::Eater(PACMAN, element-1, element); // keep moving down
			}
			if (oldelement + 1 == element) // if we came here by moving up
			{
				//cout << "DEBUG PACMAN: 2.2.2 " << endl;
				//cout << "Up" << endl;
				PacMan::Eater(PACMAN, element+1, element); // keep moving up
			}
			else if (oldelement == 9999)  // starting point for program
			{
				//cout << "DEBUG PACMAN: 2.2.3 " << endl;
				//cout << "Start" << endl;
				PacMan::Eater(PACMAN, element-1, element); // move down
				PacMan::Eater(PACMAN, element+1, element); // move up
			}
		}
		if (element == end->first && PACMAN.Elements.size() != 1) // if we are at the last element of the list but the size is larger than 1
		{
			//cout << "DEBUG PACMAN: 2.3 " << endl;
			//cout << "Down" << endl;
			PacMan::Eater(PACMAN, element-1, element); //start moving down
		}
		
		
		
	}
	//cout << "DEBUG PACMAN: 3 " << endl;
	if (PACMAN.Elements.count(element) != 0 ) // I always add the starting element and only add the new element if it exists.
	{
		//cout << "DEBUG PACMAN: 4 " << element << " " << PACMAN.Elements.size() << " " << PACMAN.Elements.count(element) << " " << oldelement << endl;
		Cluster.Elements.insert(std::pair<int,double>(element,PACMAN.Elements.at(element)));
		PACMAN.Elements.erase(element);
		//cout << "DEBUG PACMAN: 5 " << endl;
	}
	if (oldelement == 9999) 
	{
		//cout << "End" << endl;
		Cluster.SetParameters();
	}
	
};
