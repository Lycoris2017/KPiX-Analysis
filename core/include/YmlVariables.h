//-----------------------------------------------------------------------------
// File          : YmlVariables.h
// Author        : Mengqing Wu <mengqing.wu@desy.de>
// Created       : 25/10/2018
// Project       : Lycoris telescope
// Function      : Compact yml variable handler
//-----------------------------------------------------------------------------

/*
 * The Yml Data Structure: update @2018-10-23
 *  - DesyTrackerRoot as Root level
 *  - DesyTrackerRunControl/ DesyTracker/ DataWriter... as 2nd level
 * Strategy:
 *  - remove the Root level but keep 2nd level
 *  - fill all variables from 2nd level to a varHolder
 *  - for more levels under 2nd level:
 *    - map.key -> "2nd-level-name:3rd:4rd..."
 */
#ifndef __YML_VARIABLES_H__
#define __YML_VARIABLES_H__

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <map>
#include <vector>
#include <sstream> //isstringstream
#include <string>
  
/*
 * X-macro to do a mapping of the YmlLevel names
 */

#define YML_LEVEL_TABLE		\
  X(Root,       "DesyTrackerRoot")	\
  X(DataWriter, "DataWriter")		\
  X(RunControl, "DesyTrackerRunControl") \
  X(Config,     "DesyTracker") \
  X(SysLog,     "SystemLog")

#define X(a, b) a,
enum YML_LEVEL{
  YML_LEVEL_TABLE
};
#undef X

#define X(a, b)b,
// All bugs fixed by char const* instead of char*;
 // otherwise it complains string convert to char*.
 //
static char const*yml_level[] = {
  YML_LEVEL_TABLE
};
#undef X


using namespace std;

class YmlVariables{

private:
  //char* _buff;
  std::map<std::string, std::string> _vars; // !!! important -> key:value storage
  bool _debug;
  
public:

  //-- public variables
  //  const char *yml_level;

  //-- Constructor
  YmlVariables();

  //-- Deconstructor
  ~YmlVariables();

  //-- String manuplater
  string removeWhite(string str);
  
  //-- X-Macro test usage
  enum YML_LEVEL test = Root;

  //-- string processors
  //  char* c_getline(char* cin);
  
  //-- funcs
  char* fakeDataReader (const char* fn);
  bool buffParser (const char* type, char* buff); // type <- yml_level[]
  bool YmlVarReader (std::string ymlline);

  // clear the variable map
  void clear () { _vars.clear(); }
  
  void setDebug (bool debug) { _debug = debug; }

  void print() {
    puts("-- Print var map:\n");
    for ( auto& it : _vars ) cout << it.first << " : " << it.second << endl;
  }
  
  uint32_t getInt( std::string var );
  std::string getStr( std::string var);
  
};

#endif
