#!/bin/bash
analyse_files=(/scratch/data/testbeam201902/2019_03_07_16_31_32.bin /scratch/data/testbeam201902/2019_03_07_16_35_59.bin /scratch/data/testbeam201902/2019_03_07_16_37_38.bin /scratch/data/testbeam201902/2019_03_07_16_39_14.bin)
analyse_extension=.root
kpix=(17 19)


for i in ${analyse_files[*]}
do
	#python2.7 python/plot_producer.py $i$analyse_extension -n hit_diff_timed  -d 'hist e'
	python2.7 python/plot_producer.py $i$analyse_extension -n Channel_entries_k _b0 -d 'hist e' -b 0 --refuse no timed
	python2.7 python/plot_producer.py $i$analyse_extension -n Channel_entries_k _b0 -d 'hist e same' -b 0 --refuse no timed --legend S55-1 S55-2 -o S55_channel_entries.png --rebin 8
	
done



