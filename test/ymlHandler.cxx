//-----------------------------------------------------------------------------
// File          : ymlHandler.cxx
// Author        : Mengqing Wu <mengqing.wu@desy.de>
// Created       : 18/10/2018
// Update        : 24/10/2018
// Project       : read compact yml variables into a map
//-----------------------------------------------------------------------------

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <map>
#include <vector>
#include <sstream> //isstringstream
#include <string>
#include <cstring> // strchr()

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

using namespace std;

  
/*
 * X-macro to do a mapping of the YmlLevel names
 */
#define YML_LEVEL_TABLE \
  X(Root,       "DesyTrackerRoot") \
  X(DataWriter, "DataWriter") \
  X(RunControl, "DesyTrackerRunControl") \
  X(Config,     "DesyTracker")

#define X(a, b) a,
enum YML_LEVEL{
  YML_LEVEL_TABLE
};
#undef X

#define X(a, b)b,
char *yml_level[] = {
  YML_LEVEL_TABLE
};
#undef X


bool str_isInside(const std::string & str, char c){ return str.find(c) != std::string::npos;}


class ymlHandler{

private:
  //char* _buff;
  std::map<std::string, std::string> _vars; // !!! important -> key:value storage
  bool _debug;
  
public:
  //-- Constructor
  ymlHandler();

  //-- Deconstructor
  ~ymlHandler();

  //-- X-Macro test usage
  enum YML_LEVEL test = Root;

  //-- string processors
  //  char* c_getline(char* cin);
  
  //-- funcs
  char* fakeDataReader (const char* fn);
  bool buffParser (char* type, char* buff); // type <- yml_level[]
  bool YmlVarReader (std::string ymlline);

  void setDebug (bool debug) { _debug = debug; }
  void print() {
    puts("-- Print var map:\n");
    for ( auto& it : _vars ) cout << it.first << " : " << it.second << endl;
  }
  
  uint32_t getInt( std::string var );
  std::string getStr( std::string var);
  
};


/// --- application main func --- ///

int main(int argc, const char* argv[])
{
 if (argc<2){
    cout << "Usage: ymlHandler config.yml"<<endl;
    return 1;
  }
 
 ymlHandler yhand;
 auto data   = yhand.fakeDataReader( argv[1] );

 // cout << "main :: debug - data" << endl;
 // cout << data << endl;
 
 yhand.buffParser( yml_level[RunControl], data);
 yhand.print(); 
 return (0);
  
}

/// ------ class functions part: -------- ///

ymlHandler::ymlHandler(){

  if (_debug) printf("X-Macro: yml level test = %s\n", yml_level[test]);
  _vars.clear();
  _debug=false;
}

ymlHandler::~ymlHandler(){}  // empty template

uint32_t ymlHandler::getInt ( std::string var ) {
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

std::string ymlHandler::getStr( std::string var ) {

  auto iter = _vars.find( var );
  if ( iter == _vars.end() )
    return ("");
  
  if (_debug) printf(" Find your key : value -> %s : %s", iter->first, iter->second);
  return iter->second;
   return ("");
}



bool ymlHandler::YmlVarReader( std::string ymlline ){
  
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

      if ( tkKey == yml_level[Root] )
	cout << " Warning: not expect Root here   -> " << tkKey << endl;
      else if ( tkKey == yml_level[DataWriter] ||
		tkKey == yml_level[Config]     ||
		tkKey == yml_level[RunControl]
		)
	cout << " Warning: not expect Mother here -> " << tkKey << endl;
      else {// child node found!
	if ( _debug )cout << "[debug] child node found -> " << tkKey << endl;
	if ( !variable.empty() )   variable += ":";
	variable += tkKey;
      }
      
    }
    
  }

  if (_debug){
    cout << "map.first  -> " << variable << endl;;
    cout << "map.second -> " << value << endl;
  }
  _vars[variable] = value;

  return true;
}


char* ymlHandler::fakeDataReader(const char* fn) {

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

  // set the length you want to readout
  long file_length = file.tellg();
  file.seekg (0, std::ios_base::beg);
  file.clear ();

  datab = new char[file_length];
  file.read (datab, file_length);

  std:: cout << "fakeDataReader get char* : \n"<< datab << std::endl ;
  
  
  return datab;
}

bool ymlHandler::buffParser( char* type, char* buff ){

  /* funcs:
   * - loop over all the lines from the char buff[];
   * - according to " mother node ? type ", parse to VarReader or not.
    */

  std::string buff_str = buff;

  istringstream iss(buff_str);
  std::string line = "";

  if (_debug) cout << " debug: you want type -> "<< type << endl;
  
  while ( std::getline(iss, line, '\n') ){
    if ( line.empty() ) continue;

    if (_debug) cout << " Line -> " << line << endl;

    long root_pos = line.find('.');
    if ( line.substr(0, root_pos) != yml_level[Root] ){ // move to next yml level!
      if (_debug) cout << " Wrong Root -> "<< line.substr(0, root_pos);
      continue; // wrong root!
    }
    else {
      line = line.substr( root_pos+1, line.size() );
      if (_debug) cout << " Root stripped -> "<< line << endl;
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
