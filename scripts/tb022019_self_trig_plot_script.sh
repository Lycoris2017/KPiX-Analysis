#!/bin/bash
analyse_files=(/scratch/data/testbeam201902/2019_02_15_13_57_41.bin /scratch/data/testbeam201902/2019_02_15_14_12_37.bin /scratch/data/testbeam201902/2019_02_15_14_23_30.bin /scratch/data/testbeam201902/2019_02_13_11_40_07.bin /scratch/data/testbeam201902/2019_02_13_11_43_15.bin /scratch/data/testbeam201902/2019_02_13_11_46_28.bin /scratch/data/testbeam201902/2019_02_13_11_49_38.bin /scratch/data/testbeam201902/2019_02_13_11_36_31.bin /scratch/data/testbeam201902/2019_02_13_11_22_33.bin)
analyse_extension=.root
kpix=(17 19 28 30)
horiz_scan_files=(/scratch/data/testbeam201902/2019_02_13_11_49_38.bin /scratch/data/testbeam201902/2019_02_13_11_46_28.bin /scratch/data/testbeam201902/2019_02_13_11_36_31.bin  /scratch/data/testbeam201902/2019_02_13_11_40_07.bin /scratch/data/testbeam201902/2019_02_13_11_43_15.bin)
legend=(430 440 450 460 470)


for i in ${analyse_files[*]}
do
	#python python/plot_producer.py $i$analyse_extension -n hit_diff_timed  -d 'hist e'
	python python/plot_producer.py $i$analyse_extension -n right_strip -d 'hist e' -b 0 --refuse 28 30
	python python/plot_producer.py $i$analyse_extension -n timed_right -d 'hist e same' -b 0 --refuse 28 30 --legend S59-2 S43-2 -o sensor_comparison_timed_right_strip_entries.png --rebin 8
done

#echo $same
#echo $dates

same=""
legendname=""

for i in ${horiz_scan_files[*]}
do
	same="$same $i$analyse_extension"
done
for i in ${legend[*]}
do
	legendname="$legendname $i"
done

#echo $same
python python/plot_producer.py $same -n hit_diff_timed  -d 'hist e same' -o horizontal_scan.png --legend $legendname --xrange 0 40

