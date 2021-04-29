//-----------------------------------------------------------------------------
// File          : KpixSample.cpp
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 05/29/2012
// Project       : Kpix DAQ
//-----------------------------------------------------------------------------
// Description :
// Sample Container
//-----------------------------------------------------------------------------
// Copyright (c) 2012 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 05/29/2012: created
//-----------------------------------------------------------------------------
#include <string.h>
#include <iostream>
#include "KpixSample.h"
using namespace std;

// Constructor for static pointer
KpixSample::KpixSample () {
   data_        = ldata_;
   eventNumber_ = 0;
}

// Constructor with copy
KpixSample::KpixSample ( uint *data, uint eventNumber ) {
   data_        = ldata_;
   eventNumber_ = eventNumber;
   memcpy(ldata_,data,8);
}

//! DeConstructor
KpixSample::~KpixSample ( ) { }

// Set data pointer.
void KpixSample::setData ( uint *data, uint eventNumber ) {
   data_        = data;
   eventNumber_ = eventNumber;
}

// Get sample event number.
uint KpixSample::getEventNum() {
   return(eventNumber_);
}

// Get KPIX address from sample.
uint KpixSample::getKpixAddress() {
   uint ret;
   ret = (data_[0] >> 16) & 0xFFF;
   return(ret);
}

// Get KPIX channel.
uint KpixSample::getKpixChannel() {
   uint ret;
   ret = data_[0] & 0x3FF;
   return(ret);
}

// Get KPIX bucket.
uint KpixSample::getKpixBucket() {
   uint ret;
   ret = (data_[0] >> 10) & 0x3;
   return(ret);
}

// Get sample range.
uint KpixSample::getSampleRange() {
   uint ret;
   ret = (data_[0] >> 13) & 0x1;
   return(ret);
}

// Get sample time
uint KpixSample::getSampleTime() {
   uint ret;
   ret = (data_[1] >> 16) & 0x1FFF;
   return(ret);
}

// Get sample value.
uint KpixSample::getSampleValue() {
   uint ret;
   ret = data_[1] & 0x1FFF;
   return(ret);
}

// Get empty flag.
uint KpixSample::getEmpty() {
   uint ret;
   ret = (data_[0] >> 15) & 0x1;
   return(ret);
}

// Get badCount flag.
uint KpixSample::getBadCount() {
   uint ret;
   ret = (data_[0] >> 14) & 0x1;
   return(ret);
}

// Get trigger type flag.
uint KpixSample::getTrigType() {
   uint ret;
   ret = (data_[0] >> 12) & 0x1;
   return(ret);
}

// Get sample type
KpixSample::SampleType KpixSample::getSampleType() {
   SampleType ret;
   ret = (SampleType)((data_[0] >> 28) & 0xf);
   return(ret);
}


// ExtTs runtime -- last 32 bits of the FPGA clock counter
uint KpixSample::getSampleRuntime32() {
  // 63-32 bits
  uint ret;
  ret = data_[1];
  return(ret);
}

uint64_t KpixSample::getSampleRuntime64(uint64_t frameruntime, uint &overflow){
    uint64_t frameRuntime32l = (frameruntime & 0x00000000FFFFFFFF);
    uint64_t frameRuntime32h = (frameruntime & 0xFFFFFFFF00000000);
//    std::cout << overflow << endl;
//    std::cout << frameRuntime32h << endl;
    if (frameRuntime32l > data_[1]){
        overflow++;
        frameRuntime32h = (frameRuntime32h>>32)+overflow;
        frameRuntime32h = (frameRuntime32h<<32);
    }
//	uint64_t ret = (frameruntime & 0xFFFFFFFF00000000) | data_[1];
	uint64_t ret = (frameRuntime32h & 0xFFFFFFFF00000000) | data_[1];
//    std::cout << (frameRuntime32h & 0xFFFFFFFF00000000) << endl;
//    std::cout << frameRuntime32h << endl;
    return(ret);
}


// ExtTs bunchcount -- TimeStamp signal happened inside the 0-8191 BunchClkCount
uint KpixSample::getBunchCount() {
  // 15-3 bits
  uint ret;
  ret = (data_[0]>>3) & 0xFFF;
  return(ret);
}

// ExtTs bunchcount -- TimeStamp signal happened inside the 0-8191 BunchClkCount
uint KpixSample::getSubCount() {
  // 2-0 bits
  uint ret;
  ret = (data_[0] & 0x3);
  return(ret);
}


