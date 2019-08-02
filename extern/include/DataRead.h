//-----------------------------------------------------------------------------
// File          : DataRead.h
// Author        : Mengqing Wu <mengqing.wu@desy.de>
// Created       : 25/10/2018
// Project       : Lycoris Telescope DAQ
//-----------------------------------------------------------------------------
// Description :
// Read data & configuration from disk, originally from Ryan Herbst @SLAC
//-----------------------------------------------------------------------------
// Modification history :
// 21/08/2018 : scrap out of old KPiX DAQ package;
// xx/09/2018 : BZLIB free, shared memory free;
// 25/10/2018 : has Yml string from bin data to readable holder
//-----------------------------------------------------------------------------

#ifndef __DATA_READ_H__
#define __DATA_READ_H__

#include <string>
#include <map>
#include <stdint.h>

#include "Data.h"
#include "XmlVariables.h"

#include "YmlVariables.h"

/* #ifdef USE_BZLIB */
/* #include <bzlib.h> */
/* #endif */

using namespace std;

#ifdef __CINT__
#define uint32_t unsigned int
#endif

// Define variable holder
typedef map<string,string> VariableHolder;

//! Class to contain generic register data.
class DataRead {

	uint32_t rdAddr_;
	uint32_t rdCount_;
	
	// File descriptor
	int32_t fd_;
	
	// File size
	off_t size_;

	// Compression options
	bool     bzEnable_;
	BZFILE * bzFile_;
	
	// Process xml
	void xmlParse ( uint32_t size, char *data );

	// Process yml
	void ymlParse ( uint32_t size, char *data );
	
	// Variables
	XmlVariables status_;
	XmlVariables config_;
	XmlVariables start_;
	XmlVariables stop_;
	XmlVariables time_;

	YmlVariables yaml_; 
		
	// Start/Stop flags
	bool sawRunStart_;
	bool sawRunStop_;
	bool sawRunTime_;
	
 public:

	bool     debug_;

	//! Constructor
	DataRead ( );
	
	//! Deconstructor
	~DataRead ( );
	
	//! Open File
	/*! 
	 * \param file Filename
	 */
	bool open ( string file, bool compressed = false );
	
	//! Open Shared Memory
	/*! 
	 * \param system System name
	 * \param id ID to identify your process
	 */
	//void openShared ( string system, uint32_t id, int32_t uid=-1 );
	
	//! Close File
	void close ( );
	
	//! Return file size in bytes
	off_t size ( );
	
	//! Return file position in bytes
	off_t pos ( );
	
	//! Get next data record
	/*! 
	 * Returns true on success
	 * \param data Data object to store data
	 */
	bool next ( Data *data );
	
	//! Get next data record & create new data object
	/*! 
	 * Returns NULL on failure
	 */
	Data *next ( );

	string m_str_yamlst = yml_level[RunControl];
	string m_str_yamlcfg = yml_level[Config];
	//! Get a config value
	/*! 
	 * \param var Config variable name
	 */
	string getConfig ( string var );
	string getYmlConfig ( string var );

	//! Get a config value as integer
	/*! 
	 * \param var Config variable name
	 */
	uint32_t getConfigInt ( string var );
	uint32_t getYmlConfigInt ( string var );

	//! Get a status value
	/*! 
	 * \param var Status variable name
	 */
	string getStatus ( string var );
	string getYmlStatus ( string var );
	
	//! Get a status value as integer
	/*! 
	 * \param var Status variable name
	 */
	uint32_t getStatusInt ( string var );
	uint32_t getYmlStatusInt ( string var );
	
	//! Dump config
	void dumpConfig ( );
	
	//! Dump status
	void dumpStatus ( );
	
	//! Get config as XML
	string getConfigXml ( );
	
	//! Dump status
	string getStatusXml ( );
	
	//! Get a run start value
	/*! 
	 * \param var Variable name
	 */
	string getRunStart ( string var );
	
	//! Get a run stop value
	/*! 
	 * \param var Variable name
	 */
	string getRunStop ( string var );
	
	//! Get a run time value
	/*! 
	 * \param var Variable name
	 */
	string getRunTime ( string var );
	
	//! Dump start
	void dumpRunStart ( );
	
	//! Dump stop
	void dumpRunStop ( );
	
	//! Dump time
	void dumpRunTime ( );
	
	//! Return true if we saw start marker, self clearing
	bool  sawRunStart ( );
	
	//! Return true if we saw stop marker, self clearing
	bool  sawRunStop ( );
	
	//! Return true if we saw time marker, self clearing
	bool  sawRunTime ( );
};
#endif
