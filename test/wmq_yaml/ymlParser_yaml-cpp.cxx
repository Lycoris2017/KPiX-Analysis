//-----------------------------------------------------------------------------
// File          : ymlParser.cxx
// Author        : Mengqing Wu <mengqing.wu@desy.de>
// Created       : 18/10/2018
// Project       : an Yml Parser, based on yaml-cpp lib
//-----------------------------------------------------------------------------

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string.h>

#include "yaml-cpp/yaml.h"

int main()
{
  //std::ifstream fin("Test1Ch.yml");
  //YAML::Parser parser(fin);
  
  YAML::Node config = YAML::LoadFile("Test1Ch.yml");
  
  std::cout<< "Size of config : "
	   << config.size() << std::endl;
  
  switch (config.Type()) {
  case YAML::NodeType::Null:
    puts ("Null\n"); break;
  case YAML::NodeType::Scalar: // ...
    puts ("Scalar\n"); break;
  case YAML::NodeType::Sequence: // ...
    puts ("Sequence\n"); break;
  case YAML::NodeType::Map: // ...
    puts ("map\n"); break;
  default:
    printf("Got node type %d\n", config.Type());
  }
  
  return 0;
}
