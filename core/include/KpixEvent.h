//-----------------------------------------------------------------------------
// File          : KpixEvent.h
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 05/29/2012
// Project       : KPIX DAQ Software
//-----------------------------------------------------------------------------
// Description :
// Event Data consists of the following: Z[xx:xx] = Zeros
//    Header = 8 x 32-bits
//       Header[0] = EventNumber[31:0]
//       Header[1] = Timestamp[31:00] --> RunTime - LSB, little endian @MQ
//       Header[2] = Zeros[31:0] --> RunTime - HSB, little endian @MQ
//       Header[3] = Zeros[31:0]
//       Header[4] = Zeros[31:0]
//       Header[5] = Zeros[31:0]
//       Header[6] = Zeros[31:0]
//       Header[7] = Zeros[31:0]
//
//    Samples = 2 x 32-bits
//
//    Tail = 1 x 16-bits
//       Tail[0] = Zeros
//-----------------------------------------------------------------------------
// Copyright (c) 2012 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 05/29/2012: created
// 11/06/2019: modified -- Mengqing
//----------------------------------------------------------------------------
#ifndef __KPIX_EVENT_H__
#define __KPIX_EVENT_H__

#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include "KpixSample.h"
#include <Data.h>
using namespace std;

//! Tracker Event Container Class
class KpixEvent : public Data {

      // Frame Constants
      static const uint headSize_   = 8;
      static const uint tailSize_   = 1;
      static const uint sampleSize_ = 2;

      // Internal sample contrainer
      KpixSample sample_;

   public:

      //! Constructor
      KpixEvent ();

      //! Deconstructor
      ~KpixEvent ();

      //! Get event number 
      uint eventNumber ( );

      //! Get timestamp
      uint timestamp ( );

      //! Get frame runtime
      uint64_t runtime ();

      //! Get sample count
      uint count ( );

      //! Get sample at index
      /*!
       * Returns pointer to static sample object without memory allocation.
       * Contents of returned object will change next time sample() is called.
       * \param index Sample index. 0 - count()-1.
      */
      KpixSample *sample (uint index);

      //! Get sample at index
      /*!
       * Returns pointer to copy of sample object. A newly allocated sample object
       * is created and must be deleted after use.
       * \param index Sample index. 0 - count()-1.
      */
      KpixSample *sampleCopy (uint index);

};

#endif
