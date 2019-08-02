#!/bin/bash
large_file_70v=/scratch/data/testbeam201902/2019_02_13_15_37_39.bin
bias_scan_files=(/scratch/data/testbeam201902/2019_02_13_13_46_28.bin /scratch/data/testbeam201902/2019_02_13_13_50_24.bin /scratch/data/testbeam201902/2019_02_13_13_57_50.bin /scratch/data/testbeam201902/2019_02_13_14_01_29.bin /scratch/data/testbeam201902/2019_02_13_14_09_01.bin /scratch/data/testbeam201902/2019_02_13_14_13_53.bin /scratch/data/testbeam201902/2019_02_13_14_26_12.bin /scratch/data/testbeam201902/2019_02_13_14_30_39.bin /scratch/data/testbeam201902/2019_02_13_14_40_41.bin /scratch/data/testbeam201902/2019_02_13_14_45_36.bin /scratch/data/testbeam201902/2019_02_13_14_55_07.bin /scratch/data/testbeam201902/2019_02_13_14_57_20.bin /scratch/data/testbeam201902/2019_02_13_15_15_52.bin /scratch/data/testbeam201902/2019_02_13_15_26_42.bin /scratch/data/testbeam201902/2019_02_13_15_31_42.bin) 
calibration0213=/scratch/data/testbeam201902/2019_02_13_12_21_15.bin.newCalib.root

pedestal_ext=".pedestal.root"

for i in ${bias_scan_files[*]}
do
	./bin/pedestal $i $calibration0213
	./bin/analysisExternal $i $calibration0213 $i$pedestal_ext
	#./bin/analysisExternal $i $calibration0213 $pedestal0215$pedestal_ext 
done

