#!/bin/bash
files0215=(/scratch/data/testbeam201902/2019_02_15_11_22_33.bin /scratch/data/testbeam201902/2019_02_15_11_26_40.bin /scratch/data/testbeam201902/2019_02_15_11_30_22.bin /scratch/data/testbeam201902/2019_02_15_11_32_50.bin)
files0213=(/scratch/data/testbeam201902/2019_02_13_13_42_33.bin /scratch/data/testbeam201902/2019_02_13_13_46_28.bin  /scratch/data/testbeam201902/2019_02_13_15_37_39.bin /scratch/data/testbeam201902/2019_02_13_13_39_57.bin /scratch/data/testbeam201902/2019_02_13_13_36_34.bin /scratch/data/testbeam201902/2019_02_13_13_33_12.bin /scratch/data/testbeam201902/2019_02_13_13_28_11.bin /scratch/data/testbeam201902/2019_02_13_13_24_39.bin /scratch/data/testbeam201902/2019_02_13_13_19_57.bin /scratch/data/testbeam201902/2019_02_13_12_14_12.bin /scratch/data/testbeam201902/2019_02_13_12_10_35.bin /scratch/data/testbeam201902/2019_02_13_12_00_47.bin /scratch/data/testbeam201902/2019_02_13_11_56_38.bin)
fileshigh0215=(/scratch/data/testbeam201902/2019_02_15_11_30_22.bin /scratch/data/testbeam201902/2019_02_15_11_24_13.bin)

calibration0215=/scratch/data/testbeam201902/2019_02_15_10_03_05.bin.newCalib.root
calibrationhigh0215=/scratch/data/testbeam201902/2019_02_15_14_41_39.bin.newCalib.root
calibration0213=/scratch/data/testbeam201902/2019_02_13_12_21_15.bin.newCalib.root
pedestal0215=/scratch/data/testbeam201902/2019_02_15_11_22_33.bin
pedestalhigh0215=/scratch/data/testbeam201902/2019_02_15_11_24_13.bin
pedestal0213=/scratch/data/testbeam201902/2019_02_13_13_42_33.bin


pedestal_ext=".pedestal.root"

for i in ${files0215[*]}
do
	./bin/pedestal $i $calibration0215
	./bin/analysisExternal $i $calibration0215 $i$pedestal_ext
	./bin/analysisExternal $i $calibration0215 $pedestal0215$pedestal_ext 
done
for i in ${files0213[*]}
do
	./bin/pedestal $i $calibration0213
	./bin/analysisExternal $i $calibration0213 $i$pedestal_ext
	./bin/analysisExternal $i $calibration0213 $pedestal0213$pedestal_ext 
done
for i in ${fileshigh0215[*]}
do 
	./bin/pedestal $i $calibrationhigh0215
	./bin/analysisExternal $i $calibrationhigh0215 $i$pedestal_ext
	./bin/analysisExternal $i $calibrationhigh0215 $pedestalhigh0215$pedestal_ext
done

