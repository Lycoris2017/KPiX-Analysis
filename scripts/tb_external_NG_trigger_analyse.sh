#!/bin/bash

extension="_external.root"

#unused files /scratch/data/tracker_test/2018_08_17_15_31_33.bin /scratch/data/tracker_test/eudaq/testbeam_201808/2018_08_17_20_30_20.bin /scratch/data/tracker_test/eudaq/testbeam_201808/2018_08_17_20_01_21.bin  /scratch/data/tracker_test/2018_08_17_15_31_33.bin
data_files_normal_gain_320ns=(/scratch/data/tracker_test/2018_08_19_21_34_17.bin /scratch/data/tracker_test/2018_08_17_15_31_33.bin)
pedestal_normal_gain=/scratch/data/tracker_test/2018_08_19_21_34_17.bin.pedestal.root
datapedestal_normal_gain=/scratch/data/tracker_test/2018_08_17_15_31_33.bin.pedestal.root
calibration_file_normal_gain_320ns=/scratch/data/tracker_test/2018_08_19_20_16_08.bin.newCalib.root
for i in ${data_files_normal_gain_320ns[*]}
do
	#./bin/analysisExternal $i $calibration_file_normal_gain_320ns $pedestal_normal_gain
	./bin/analysisExternal $i $calibration_file_normal_gain_320ns $datapedestal_normal_gain
done



#for i in ${data_files_normal_gain_320ns[*]}
#do 
	#echo "$i$extension"
#done
#data_files_high_gain_320ns=(/scratch/data/tracker_test/eudaq//testbeam_201808/2018_08_18_16_31_51.bin /scratch/data/tracker_test/2018_08_19_21_54_13.bin)
#pedestal_high_gain=/scratch/data/tracker_test/2018_08_19_21_54_13.bin.root
#calibration_file_high_gain_320ns=/scratch/data/tracker_test/2018_09_05_11_06_13.bin.newCalib.root
#for i in ${data_files_high_gain_320ns[*]}
#do
	#./bin/analysisExternal $i $calibration_file_high_gain_320ns $pedestal_high_gain
#done
