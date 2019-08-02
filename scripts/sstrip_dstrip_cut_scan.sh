#!/bin/bash

files=(/scratch/data/tracker_test/2018_08_17_15_31_33.bin_2018_08_17_15_31_33.pedestal.external.root)


for i in ${files[*]}
do
	python python/plot_producer.py $i -n cuts_singlestrip -b 0 -k 26 -d "h e" --xrange 0 20 -o "cuts2.1_singlestrip_0_20.png"
	python python/plot_producer.py $i -n cuts_doublestrip -b 0 -k 26 -d "h e" --xrange 0 20 -o "cuts0.85_doublestrip_0_20.png"
done
