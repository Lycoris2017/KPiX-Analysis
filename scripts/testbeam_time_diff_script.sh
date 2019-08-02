#!/bin/bash

#unused files /scratch/data/tracker_test/eudaq//testbeam_201808/2018_08_18_20_21_17.bin
data_files=( /scratch/data/tracker_test/eudaq//testbeam_201808/2018_08_18_20_21_17.bin /scratch/data/tracker_test/eudaq//testbeam_201808/2018_08_18_20_06_16.bin /scratch/data/tracker_test/eudaq//testbeam_201808/2018_08_18_16_48_58.bin /scratch/data/tracker_test/eudaq//testbeam_201808/2018_08_18_16_31_51.bin /scratch/data/tracker_test/eudaq/testbeam_201808/2018_08_20_00_51_01.bin /scratch/data/tracker_test/2018_08_17_10_32_06.bin) #(/scratch/data/tracker_test/2018_03_28_16_13_28.bin /scratch/data/tracker_test/2018_03_28_16_18_31.bin)
analysis_extension=.root
calibration_file=/scratch/data/tracker_test/2018_08_19_20_16_08.bin.newCalib.root
for i in ${data_files[*]}
do
	./bin/analysis $i /scratch/data/tracker_test/2018_08_19_20_16_08.bin.newCalib.root
	#python python/plot_producer.py $i$analysis_extension -d 'same he' -n ext_time_diff --refuse con -o time_diff_zoom_ylog.png --ylog True --yrange 5E-3 4 --xrange -20 20 --legend Total Tracker ECAL
	#python python/plot_producer.py $i$analysis_extension -d 'same he' -n ext_time_diff --refuse con -o time_diff_zoom.png --xrange -20 20 --legend Total Tracker ECAL
	#python python/plot_producer.py $i$analysis_extension -d 'same he' -n ext_time_diff --refuse con -o time_diff.png --xrange -150 150 --legend Total Tracker ECAL
	#python python/plot_producer.py $i$analysis_extension -n timestamp_kpix_k_ --refuse monster -b 0 4 -d 'hist e'
	#python python/plot_producer.py $i$analysis_extension -n Left_Strip_entries_k_26 -b 0 4  -d 'hist e'
done





