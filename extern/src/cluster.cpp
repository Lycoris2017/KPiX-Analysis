#include <iostream>
#include <string>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iterator>

#include "cluster.h"
using namespace std;

void cluster::construct(vector<int> hit_position, vector<double> hit_charge)
{
	double weighted_position = 0;
	double charge_sum = 0;
	
	
	int vec_mxch = distance(hit_charge.begin(), max_element( hit_charge.begin(), hit_charge.end() )); // position of the highest hit charge element in the vector
	double maxcharge = hit_charge.at(vec_mxch);
	double maxcharge_position = hit_position.at(vec_mxch);
	
	clustering(hit_position, hit_charge, maxcharge_position, 1);
	clustering(hit_position, hit_charge, maxcharge_position, -1);
	clustercharge -= maxcharge;
	
	//for (int i = 0; i < hit_position.size(); i++)
	//{
		//if (hit_position.at(i) < maxcharge_position+2 && hit_position.at(i) > maxcharge_position-2)
		//{
			//weighted_position += hit_charge.at(i)*hit_position.at(i);
			//charge_sum += hit_charge.at(i);
		//}
	//}
	//position = weighted_position/charge_sum;
	//charge = cluster;
	
	position = 1;
	charge = 1;
}


void cluster::clustering(vector<int> hit_position, vector<double> hit_charge, int start, int step)
{
	if (hit_position.at(start+step) == hit_position.at(start)+step)
	{
		clustercharge += hit_charge.at(start);
		clustering(hit_position, hit_charge, hit_position.at(start)+step, step);
	}
	else
	{
		clustercharge += hit_charge.at(start);
	}
}
	
	


