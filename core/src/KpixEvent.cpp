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
//       Header[1] = Timestamp[31:00] --> runtime lsb @MQ
//       Header[2] = Zeros[31:0] --> runtime hsb @MQ
//       Header[3] = Zeros[31:0]
//       Header[4] = Zeros[31:0]
//       Header[5] = Zeros[31:0]
//       Header[6] = Zeros[31:0]
//       Header[7] = Zeros[31:0]
//
//    Samples = 2 x 32-bits
//
//    Tail = 1 x 32-bits
//       Tail[0] = Zeros
//       Tail[1] = Zeros
//-----------------------------------------------------------------------------
// Copyright (c) 2012 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 05/29/2012: created
//----------------------------------------------------------------------------
// Modified by Mengqing Wu for Lycoris Project @ 28/09/2018
// Modified by MQ for Lycoris @ 11/06/2019
// -- 
//----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "KpixEvent.h"
#include <stdint.h>
using namespace std;

// Constructor
KpixEvent::KpixEvent () : Data() { }

// Deconstructor
KpixEvent::~KpixEvent () { }

// Get event number
uint KpixEvent::eventNumber ( ) {
   return(data_[0]);
}

// Get timestamp
uint KpixEvent::timestamp ( ) {
   return(data_[1]);
}


// Get 64bit runtime
uint64_t KpixEvent::runtime ( ) {
    uint64_t lsb = data_[1];
    uint64_t hsb = data_[2];

    uint64_t ret =  hsb<<32 | lsb;
	
    // printf("Counter: %lu %lu; ret : %lu\n", hsb,lsb,ret); // debug
	
	return (ret);
}

// Get sample count
uint KpixEvent::count ( ) {
   uint rem = 0;
   if ( size_ <= (headSize_ + tailSize_)) return(0);

   rem = (size_-(headSize_ + tailSize_));

   if ( (rem % sampleSize_) != 0 ) return(0);

   return(rem/sampleSize_);
}

// Get sample at index
KpixSample *KpixEvent::sample (uint index) {
   if ( index >= count() ) return(NULL);
   else {
      sample_.setData(&(data_[headSize_+(index*sampleSize_)]),eventNumber());
      return(&sample_);
   }
}

// Get sample at index
KpixSample *KpixEvent::sampleCopy (uint index) {
   KpixSample *tmp;

   if ( index >= count() ) return(NULL);
   else {
      tmp = new KpixSample (&(data_[headSize_+(index*sampleSize_)]),eventNumber());
      return(tmp);
   }
}

