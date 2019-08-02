//-----------------------------------------------------------------------------
// File          : ymlEmitter.cxx
// Author        : Mengqing Wu <mengqing.wu@desy.de>
// Created       : 18/10/2018
// Project       : an Yml emitter, to interpret the compact yml format to flat yml, based on yaml-cpp lib.
//-----------------------------------------------------------------------------

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string.h>

#include "yaml-cpp/yaml.h"

int main()
{

  YAML::Emitter out;
  // out<< "Hello, World!";

  // std::cout << "Here is the output YAML:\n"
  // 	    << out.c_str() << std::endl;


  std::vector <int> squares;
  squares.push_back(1);
  squares.push_back(4);
  squares.push_back(9);
  squares.push_back(16);
  
  std::map <std::string, int> ages;
  ages["Daniel"] = 26;
  ages["Jesse"] = 24;
  
  //  YAML::Emitter out;
  out << YAML::BeginSeq;
  out << YAML::Flow << squares;
  out << ages;
  out << YAML::EndSeq;

  std::cout << "Here is the output YAML:\n"
   	    << out.c_str() << std::endl;;
  
  return 0;
}


