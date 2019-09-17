//-----------------------------------------------------------------------------
// File          : YmlVariables.cpp
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

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <map>
#include <vector>
#include <sstream> //isstringstream
#include <string>

#include "YmlVariables.h"

using namespace std;
 
// #define X(a, b)b,
//   /* All bugs fixed by char const* instead of char*;
//    * otherwise it complains string convert to char*.
//    */
//   char const* yml_level[] = {
//     YML_LEVEL_TABLE
//   };
// #undef X

/// ------ class functions part: -------- ///

YmlVariables::YmlVariables(){

  //yml_level = _yml_level;
  if (_debug) printf("X-Macro: yml level test = %s\n", yml_level[test]);
  _vars.clear();
  _debug=false;
  
}

YmlVariables::~YmlVariables(){}  // empty template

// Remove whitespace and newlines -- copy from Ryan's XmlVariable.cpp
string YmlVariables::removeWhite ( string str ) {
   string temp;
   uint32_t   i;

   temp = "";

   for (i=0; i < str.length(); i++) {
      if ( str[i] != ' ' && str[i] != '\n' ) 
         temp += str[i];
   }
   return(temp);
}


uint32_t YmlVariables::getInt ( std::string var ) {
  uint32_t     ret;
  string       value;
  const char   *sptr;
  char         *eptr;
  
  auto iter = _vars.find( var );

  if ( iter == _vars.end() ) return (0);

  value = iter->second;
  sptr = value.c_str();
  ret = (uint32_t)strtoul(sptr, &eptr, 0);
  if ( *eptr != '\0' || eptr == sptr ) ret = 0;

  return (ret); 
}

std::string YmlVariables::getStr( std::string var ) {

  auto iter = _vars.find( var );
  if ( iter == _vars.end() ) return ("");

  if (_debug)
    printf(" Find your key : value -> %s : %s\n", (iter->first).c_str(), (iter->second).c_str() );
  
  return iter->second;

}



bool YmlVariables::YmlVarReader( std::string ymlline ){
  
  /* funcs:
   * - a valid char* _buff --> read "Root.Mom.C1.C2...Var:Val" to a _vars["C1:C2...:Var"]=Val;
  */
  
  if ( ymlline.empty() )
    return false;
  
  istringstream iss(ymlline);

  if (_debug)
    std::cout << " RX yml line stream : " << iss.str() << endl;
  
  // std::vector<std::string> tkKeys;
  std::string tkKey;
  std::string variable = "";
  std::string value = "";

  while (std::getline(iss, tkKey, '.')) {
    if ( tkKey.empty() ) continue;

    //-- find the "name:value", then split them:
    std::size_t found = tkKey.find(':');
    if ( found != std::string::npos ){
      if (_debug){
	cout << " variable found: " << tkKey << endl;
	cout << " name  : " << tkKey.substr( 0, found) << endl;
	cout << " value : " << tkKey.substr( found+1, tkKey.size() ) << endl;
      }
      value =  tkKey.substr( found+1, tkKey.size() );

      if ( !variable.empty() )   variable += ":";
      variable += tkKey.substr( 0, found);
      
    }
    //-- get all levels name stacked:
    else{ 
      //tkKeys.push_back(tkKey);
      /*
      if ( tkKey == yml_level[Root] )
	cout << " Warning: not expect Root here   -> " << tkKey << endl;
      else if ( tkKey == yml_level[DataWriter] ||
		tkKey == yml_level[Config]     ||
		tkKey == yml_level[RunControl]
		)
	cout << " Warning: not expect Mother here -> " << tkKey << endl;
      
      else {// child node found!
      */
	if ( _debug )cout << "[debug] child node found -> " << tkKey << endl;
	if ( !variable.empty() )   variable += ":";
	variable += tkKey;
	/*}*/
      
    }
    
  }

  if (_debug){
    cout << "map.first  -> " << variable << endl;;
    cout << "map.second -> " << value << endl;
  }

  if ( removeWhite(value) != "" )  value = removeWhite(value);
  _vars[variable] = value;

  return true;
}


char* YmlVariables::fakeDataReader(const char* fn) {

  /* funcs:
   * a fake DataRead wrapper to get a char* _buff readout as if from DataReader;
   */

  char* datab = nullptr;
   
  std::ifstream file;
  file.open(fn, std::ios_base::in|std::ios_base::ate);
  
  if (!file.is_open()) {
    cout<< "Error: in file not open" << endl;
    return nullptr;
  }

  // sethe length you want to readout
  long file_length = file.tellg();
  file.seekg (0, std::ios_base::beg);
  file.clear ();

  datab = new char[file_length];
  file.read (datab, file_length);

  std:: cout << "fakeDataReader get char* : \n"<< datab << std::endl ;
  
  
  return datab;
}

bool YmlVariables::buffParser( const char* type, char* buff ){

  /* funcs:
   * - loop over all the lines from the char buff[];
   * - according to " mother node ? type ", parse to VarReader or not.
    */

  bool stripRoot = ( strcmp(type, yml_level[Root]) != 0 ) ;
    
  std::string buff_str = buff;

  istringstream iss(buff_str);
  std::string line = "";

  if (_debug) cout << " debug: you want type -> "<< type << endl;

  while ( std::getline(iss, line, '\n') ){
    if ( line.empty() ) continue;

    if (_debug) cout << " Line -> " << line << endl;

    if ( stripRoot ){
      long root_pos = line.find('.');
      if ( line.substr(0, root_pos) != yml_level[Root] ){ // move to next yml level!
	if (_debug) cout << " Wrong Root -> "<< line.substr(0, root_pos);
	continue; // wrong root!
      }
      else {
	line = line.substr( root_pos+1, line.size() );
	if (_debug) cout << " Root stripped -> "<< line << endl;
      }
    }
    
    long mom_pos = line.find('.');
    if ( (line.substr(0, mom_pos)) == type ) {// You find MOM!
      if (_debug) 
	cout << "Find mom! -> " << line.substr(0, mom_pos) << "\n"
	     << " debug: ymlline -> " << line.substr(mom_pos+1, line.size() ) << endl;
      
      YmlVarReader( line.substr(mom_pos+1, line.size() ) );
    }
    else {
      if (_debug) cout << " Wrong Mom! -> " <<  line.substr(0, mom_pos) << endl;
      continue; // wrong mom, move to next!
    }

  }

  return true;
 
}



///-----------------------------------------------------///

void nonsense () {

  // input:
  std::string fn = "calstate.yml";
  uint  iline = 1;

  // start:
  std::ifstream ifs;
  
  ifs.open(fn, std::ifstream::in);
  if (!ifs.is_open()) {
    cout<< "Error: infile not open" << endl;
    return ;
  }
  
  ifs.seekg(0, ios::beg);

  std::string sline;
  uint nline = 0;
  
  while ( ifs.good() ) {
    nline++;
    std::getline( ifs, sline);

    std::cout << "line is : " <<  sline << std::endl;
	
    if ( nline==iline ) 
      break;
    
  }

  /* 
   * Allocate memory dynamically and copy the content of the original string. 
   * The memory will be valid until you explicitly release it using "free". 
   * Forgetting to release it results in memory leak.
   */
  
  char* buff = (char *)malloc(sline.size() + 1);
  memcpy( buff, sline.c_str(), sline.size() + 1 );

  std:: cout << "char* : "<< buff << std::endl;

  ifs.close();

}
