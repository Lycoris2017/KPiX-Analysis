#!/bin/bash
calibration=(/scratch2/data/testbeam201902/2019_02_08_14_29_17.bin
/scratch2/data/testbeam201902/2019_02_08_16_26_23.bin
/scratch2/data/testbeam201902/2019_02_11_10_43_00.bin
/scratch2/data/testbeam201902/2019_02_12_17_19_18.bin
/scratch2/data/testbeam201902/2019_02_13_12_21_15.bin
/scratch2/data/testbeam201902/2019_02_14_11_54_29.bin
/scratch2/data/testbeam201902/2019_02_15_10_03_05.bin
/scratch2/data/testbeam201902/2019_02_20_10_41_07.bin)
#tdd0=/scratch2/data/testbeam201902/2019_02_12_17_45_10.bin 
calibration_high=(/scratch2/data/testbeam201902/2019_02_15_14_41_39.bin  /scratch2/data/testbeam201902/2019_02_11_11_36_19.bin /scratch2/data/testbeam201902/2019_02_07_12_00_47.bin)
xml=/home/lycoris-dev/kpixDaq/kpix.master/xml/testbeam201902/finalXMLs/cal320_S43-2_slot19_S59-2_slot17_ECALnew1_slot28_ECALW_slot30_TDD8.xml
xml_high=/home/lycoris-dev/kpixDaq/kpix.master/xml/testbeam201902/finalXMLs/cal320_S43-2_slot19_S59-2_slot17_ECALnew1_slot28_ECALW_slot30_TDD8_HGT.xml



for i in ${calibration[*]}
do
	./bin/new_CalibrationFitter $xml $i
done
for i in ${calibration_high[*]}
do
	./bin/new_CalibrationFitter $xml_high $i
done


