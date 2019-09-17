//-----------------------------------------------------------------------------
// File          : Data.h
// Author        : Mengqing Wu <mengqing.wu@desy.de>
// Created       : 25/10/2018
// Project       : Lycoris Telescope DAQ
//-----------------------------------------------------------------------------
// Description :
// Generic data container, originally from Ryan Herbst @SLAC
//-----------------------------------------------------------------------------
// Modification history :
// 25/10/2018 :  add Yml access
//-----------------------------------------------------------------------------

#ifndef __DATA_H__
#define __DATA_H__

#include <string>
#include <stdint.h>

#ifdef USE_BZLIB
#include <bzlib.h>
#else
#include <stdio.h>
#define BZFILE FILE
#endif

using namespace std;

#ifdef __CINT__
#define uint32_t unsigned int
#endif

//! Class to contain generic register data.
class Data {

      // Allocation
      uint32_t alloc_;

   protected:

      // Data container
      uint32_t *data_;

      // Size value
      uint32_t size_;

      // Update frame state
      virtual void update();

   public:

      // Data types. 
      // Count is n*32bits for type = 0, byte count for all others
      enum DataType {
         RawData     = 0,
         XmlConfig   = 1,
         XmlStatus   = 2,
         XmlRunStart = 3,
         XmlRunStop  = 4,
         XmlRunTime  = 5,
         YmlConfig   = 6
      };

      //! Constructor
      /*! 
       * \param data Data pointer
       * \param size Data size
      */
      Data ( uint32_t *data, uint32_t size );

      //! Constructor
      Data ();

      //! Deconstructor
      virtual ~Data ();

      //! Read data from file descriptor
      /*! 
       * \param fd File descriptor
       * \param size Data size
      */
      bool read ( int32_t fd, uint32_t size );

      //! Read data from compressed file
      /*! 
       * \param bzFile Compressed file 
       * \param size Data size
      */
      bool read ( BZFILE *bzFile, uint32_t size );

      //! Copy data from buffer
      /*! 
       * \param data Data pointer
       * \param size Data size
      */
      void copy ( uint32_t *data, uint32_t size );

      //! Get pointer to data buffer
      uint32_t *data ( );

      //! Get data size
      uint32_t size ( );

};
#endif
