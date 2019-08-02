#!/bin/bash
files0215=(/scratch/data/testbeam201902/2019_02_15_13_57_41.bin /scratch/data/testbeam201902/2019_02_15_14_12_37.bin /scratch/data/testbeam201902/2019_02_15_14_23_30.bin)
files0213=(/scratch/data/testbeam201902/2019_02_13_11_40_07.bin /scratch/data/testbeam201902/2019_02_13_11_43_15.bin /scratch/data/testbeam201902/2019_02_13_11_46_28.bin /scratch/data/testbeam201902/2019_02_13_11_49_38.bin /scratch/data/testbeam201902/2019_02_13_11_36_31.bin /scratch/data/testbeam201902/2019_02_13_11_22_33.bin /scratch/data/testbeam201902/2019_02_13_11_53_08.bin /scratch/data/testbeam201902/2019_02_13_12_03_53.bin /scratch/data/testbeam201902/2019_02_13_12_07_24.bin /scratch/data/testbeam201902/2019_02_13_12_17_20.bin /scratch/data/testbeam201902/2019_02_13_13_17_37.bin /scratch/data/testbeam201902/2019_02_13_13_22_26.bin /scratch/data/testbeam201902/2019_02_13_13_26_27.bin /scratch/data/testbeam201902/2019_02_13_13_31_33.bin /scratch/data/testbeam201902/2019_02_13_13_34_57.bin /scratch/data/testbeam201902/2019_02_13_13_38_20.bin /scratch/data/testbeam201902/2019_02_13_13_44_43.bin /scratch/data/testbeam201902/2019_02_13_13_48_41.bin /scratch/data/testbeam201902/2019_02_13_13_56_14.bin /scratch/data/testbeam201902/2019_02_13_13_59_49.bin /scratch/data/testbeam201902/2019_02_13_14_07_19.bin /scratch/data/testbeam201902/2019_02_13_14_12_05.bin /scratch/data/testbeam201902/2019_02_13_14_24_28.bin /scratch/data/testbeam201902/2019_02_13_14_28_55.bin /scratch/data/testbeam201902/2019_02_13_14_38_58.bin /scratch/data/testbeam201902/2019_02_13_14_43_43.bin /scratch/data/testbeam201902/2019_02_13_14_52_44.bin /scratch/data/testbeam201902/2019_02_13_15_14_03.bin /scratch/data/testbeam201902/2019_02_13_15_22_50.bin /scratch/data/testbeam201902/2019_02_13_15_24_58.bin)
calibration0215=/scratch/data/testbeam201902/2019_02_15_10_03_05.bin.newCalib.root
calibration0213=/scratch/data/testbeam201902/2019_02_13_12_21_15.bin.newCalib.root




for i in ${files0213[*]}
do
	./bin/analysis $i $calibration0213
done

for i in ${files0215[*]}
do
	./bin/analysis $i $calibration0215
done

