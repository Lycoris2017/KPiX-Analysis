//-----------------------------------------------------------------------------
// File          : KpixSample.h
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 05/29/2012
// Project       : Kpix DAQ
//-----------------------------------------------------------------------------
// Description :
// Sample Container
//    Samples = 2 * 32-bits
//       Sample[0] = Type[3:0],KpixId[11:0],E[0],C[0],R[0],T[0],Bucket[1:0],Chan[9:0]
//       Sample[1] = Zeros[2:0],Time[12:0],Zeros[2:0],AdcValue[12:0]
//       T = Trigger bit, 1 = external trigger
//       R = Range bit, '1' = low gain
//       C = Bad count flag
//       E = Empty sample bit
//-----------------------------------------------------------------------------
// Copyright (c) 2012 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 05/29/2012: created
// 04/05/2019: Modified for new DAQ board data stream -- Mengqing
//-----------------------------------------------------------------------------
#ifndef __KPIX_SAMPLE_H__
#define __KPIX_SAMPLE_H__
#include <sys/types.h>
using namespace std;

#ifdef __CINT__
#define uint unsigned int
#endif

//! Tracker Event Container Class
class KpixSample {

   public:

      // Sample types
      enum SampleType {
         Data        = 0,
         Temperature = 1, // kpix on-chip 
         Timestamp   = 2, 
         Runtime     = 3  // kpix data sent out time
      };

   private:

      // Local data
      uint ldata_[2];

      // Data pointer
      uint *data_;

      // Event number
      uint eventNumber_;

   public:

      //! Constructor for static pointer
      KpixSample ();

      //! Constructor with copy
      KpixSample ( uint *data, uint eventNumber );

      //! DeConstructor
      ~KpixSample ( );

      //! Set data pointer.
      /*!
       * \param data Data pointer.
      */
      void setData ( uint *data, uint eventNumber );

      //! Get sample event number.
      /*!
         Method to return sample event number. This serial number associates this sample 
         with an acquisition cycle, also refered to as a bunch train in the ILC.
         \return Sample Event Number.
         \see trainNum
      */
      uint getEventNum();

      //! Get KPIX address from sample.
      /*!
         Method to return KPIX address from sample. This is the physical address of the KPIX 
         which generated this sample.  The root file which stores samples will contain a tree of 
         KPIX objects, one for each KPIX in the system. This address can be used to associate
         this sample with a given KPIX in the tree.
         \see kpixAddress
         \return KPIX Address.
      */
      uint getKpixAddress();

      //! Get KPIX channel.
      /*!
         Method to return KPIX channel number from sample.  The KPIX contains 1024 pixels organized in an 
         array of 32 rows x 32 columns. Bits [4:0] of the channel value define the row, while bits [9:5]
         define the column.
         \see kpixChannel
         \return Kpix Channel Number, masked to ensure it lies between 0 and 1023.
      */
      uint getKpixChannel();

      //! Get KPIX bucket.
      /*!
         Method to return KPIX bucket number from sample.  Each pixel in the KPIX device can store 
         up to 4 samples. Each sample is stored in a location refered to as a bucket. This value 
         defines which bucket, numbered 0-3, contained this sample.
         \see kpixBucket
         \return Kpix Bucket Number, masked to ensure it lies between 0 and 3.
      */
      uint getKpixBucket();

      //! Get sample range.
      /*!
         The range flag determines if the sample is normal gain or low gain.
         \see sampleRange
         \return Sample Range.
               - 0 = Normal
               - 1 = Low Gain
               .
      */
      uint getSampleRange();

      //! Get sample time
      /*!
         This sample time value contains a timestamp identifying the relative
         time of the sample within the bunch train. Each bunch period, or crossing,
         in a bunch train has a unique time value. Each train can contain between
         1 and 8191 bunch crossing periods, each one with a period of 8 acquisition
         clock cycles.
         \see sampleTime
         \return Sample Time, masked to ensure it lies between 0 and 8191.
      */
      uint getSampleTime();

      //! Get sample value.
      /*!
         The sample value is the ADC value of the sample, between 0 and 8191. 
         \see sampleValue
         \return Sample value, masked to ensure it lies between 0 and 8191.
      */
      uint getSampleValue();

      //! Get empty flag.
      /*!
         The empty flag defines that the sample does not contain real data. This may occur if
         the run was performed in a mode forcing all sample values to be returned, regardless of
         the channel's event count.
         \see sampleRange
         \return Empty Flag
            - 0 = Normal
            - 1 = Empty Record
            .
      */
      uint getEmpty();

      //! Get badCount flag.
      /*!
         The bad count flag is set when the channel's event count field contained an invalid value.
         This may occur if the run was performed in a mode forcing all sample values to be returned, 
         regardless of the channel's event count.
         \see sampleRange
         \return Bad Count Flag.
            - 0 = Normal
            - 1 = Bad Count
            .
      */
      uint getBadCount();

      //! Get trigger type flag.
      /*!
         The trigger type flag distinguishes between samples caused by the channel's local trigger
         logic or an external source. The external source can either be the external force trigger
         input or the neighboring channel if nearest neighbor triggering is enabled.
         \see sampleRange
         \return Trigger Type Flag.
            - 0 = External
            - 1 = Local
            .
      */
      uint getTrigType();

      //! Get sample type
      /*!
         \return Sample type
      */
      SampleType getSampleType();


      /*
	New Data Time stamp structure -- Mengqing
       */
      uint getSampleRuntime32(); // runtime lsb based on 64-bit frame counter
      uint64_t getSampleRuntime64(uint64_t frameruntime); // you need to input the kpixEvent level runtime and the 32 bit overflow counter which is incremanted inside;
      uint getBunchCount(); // TimeStamp bunch clock count
      uint getSubCount(); // the 8 sub counts under 1 bunch count

      /*
	New RunTime structure -- Mengqing
      */

      
};

#endif
