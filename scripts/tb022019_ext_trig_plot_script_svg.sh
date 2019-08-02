#!/bin/bash
files0215=(/scratch/data/testbeam201902/2019_02_15_11_22_33.bin /scratch/data/testbeam201902/2019_02_15_11_26_40.bin /scratch/data/testbeam201902/2019_02_15_11_30_22.bin /scratch/data/testbeam201902/2019_02_15_11_32_50.bin)
files0213=(/scratch/data/testbeam201902/2019_02_13_13_42_33.bin /scratch/data/testbeam201902/2019_02_13_13_46_28.bin  /scratch/data/testbeam201902/2019_02_13_15_37_39.bin /scratch/data/testbeam201902/2019_02_13_13_39_57.bin /scratch/data/testbeam201902/2019_02_13_13_36_34.bin /scratch/data/testbeam201902/2019_02_13_13_33_12.bin /scratch/data/testbeam201902/2019_02_13_13_28_11.bin /scratch/data/testbeam201902/2019_02_13_13_24_39.bin /scratch/data/testbeam201902/2019_02_13_13_19_57.bin /scratch/data/testbeam201902/2019_02_13_12_14_12.bin /scratch/data/testbeam201902/2019_02_13_12_10_35.bin /scratch/data/testbeam201902/2019_02_13_12_00_47.bin /scratch/data/testbeam201902/2019_02_13_11_56_38.bin)

vertical_scan=(/scratch/data/testbeam201902/2019_02_13_13_39_57.bin /scratch/data/testbeam201902/2019_02_13_13_36_34.bin /scratch/data/testbeam201902/2019_02_13_13_33_12.bin /scratch/data/testbeam201902/2019_02_13_13_28_11.bin /scratch/data/testbeam201902/2019_02_13_13_24_39.bin /scratch/data/testbeam201902/2019_02_13_13_19_57.bin /scratch/data/testbeam201902/2019_02_13_12_14_12.bin /scratch/data/testbeam201902/2019_02_13_12_10_35.bin /scratch/data/testbeam201902/2019_02_13_12_00_47.bin /scratch/data/testbeam201902/2019_02_13_11_56_38.bin /scratch/data/testbeam201902/2019_02_13_13_46_28.bin)

pedestal0215=2019_02_15_11_22_33
pedestal0213=2019_02_13_13_42_33

extension=.pedestal.external.root

for i in ${files0215[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	python2.7 python/plot_producer.py $i"_"$substring$extension -n slopes_vs_stripr -b 0 -d 'hist e'
	python2.7 python/plot_producer.py $i"_"$substring$extension -n slopes_vs_channel -b 0 -d 'hist e'
	
	python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n hit_position_singlestrip_r_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n hit_position_singlestrip_r_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n hit_position_doublestrip_r_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n hit_position_doublestrip_r_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n hit_position_single+double_strip_r_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n hit_position_single+double_strip_r_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n fc_response_median_made _k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_median_made _k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n fc_response_cuts_singlestrip_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_cuts_singlestrip_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n fc_response_cuts_doublestrip_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_cuts_doublestrip_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n cluster_position -b 0 -d 'h e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n cluster offset -d 'h e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n cluster_charge -b 0 -d 'h e' --refuse 28 30
	
done

for i in ${files0213[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	python2.7 python/plot_producer.py $i"_"$substring$extension -n slopes_vs_stripr -b 0 -d 'hist e'
	python2.7 python/plot_producer.py $i"_"$substring$extension -n slopes_vs_channel -b 0 -d 'hist e'
	
	python2.7 python/plot_producer.py $i"_"$pedestal0213$extension -n hit_position_singlestrip_r_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n hit_position_singlestrip_r_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0213$extension -n hit_position_doublestrip_r_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n hit_position_doublestrip_r_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0213$extension -n hit_position_single+double_strip_r_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n hit_position_single+double_strip_r_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0213$extension -n fc_response_cuts_singlestrip_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_cuts_singlestrip_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0213$extension -n fc_response_cuts_doublestrip_k -b 0 -d 'hist e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_cuts_doublestrip_k -b 0 -d 'hist e' --refuse 28 30
	
	python2.7 python/plot_producer.py $i"_"$pedestal0213$extension -n cluster_position -b 0 -d 'h e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$pedestal0213$extension -n cluster offset -d 'h e' --refuse 28 30
	python2.7 python/plot_producer.py $i"_"$pedestal0213$extension -n cluster_charge -b 0 -d 'h e' --refuse 28 30
	
done

same=""

legends=(702.7 707.7 712.7 717.7 722.7 752.7 747.7 742.7 737.7 732.7 727.7)
count=0
for i in ${vertical_scan[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	same="$same $i"_"$substring$extension"
	python2.7 python/plot_producer.py $i"_"$substring$extension -n hit_position_single+double_strip_r_k17 -b 0 -d 'hist e same' --refuse 28 30 --legend ${legends[count]} -o k17_vertical_scan_single+double_hit_position_${legends[count]} 
	python2.7 python/plot_producer.py $i"_"$substring$extension -n hit_position_single+double_strip_r_k19 -b 0 -d 'hist e same' --refuse 28 30 --legend ${legends[count]} -o k19_vertical_scan_single+double_hit_position_${legends[count]}
	python2.7 python/plot_producer.py $i"_"$substring$extension -n cluster_position k17 -b0 -d 'h e same' --refuse 28 30 --legend ${legends[count]} -o k17_vertical_scan_cluster_position_${legends[count]}
	python2.7 python/plot_producer.py $i"_"$substring$extension -n cluster_position k19 -b0 -d 'h e same' --refuse 28 30 --legend ${legends[count]} -o k19_vertical_scan_cluster_position_${legends[count]} 
	((count++))
done


python2.7 python/plot_producer.py $same -n hit_position_single+double_strip_r_k17 -b 0 -d 'hist e same' --refuse 28 30 --legend 702.7 707.7 712.7 717.7 722.7 752.7 747.7 742.7 737.7 732.7 727.7 -o k17_vertical_scan_single+double_hit_position
python2.7 python/plot_producer.py $same -n hit_position_single+double_strip_r_k19 -b 0 -d 'hist e same' --refuse 28 30 --legend 702.7 707.7 712.7 717.7 722.7 752.7 747.7 742.7 737.7 732.7 727.7 -o k19_vertical_scan_single+double_hit_position
python2.7 python/plot_producer.py $same -n cluster_position k17 -b0  -d 'hist e same' --refuse 28 30 --legend 702.7 707.7 712.7 717.7 722.7 752.7 747.7 742.7 737.7 732.7 727.7 -o k17_vertical_scan_cluster_position
python2.7 python/plot_producer.py $same -n cluster_position k19 -b0  -d 'hist e same' --refuse 28 30 --legend 702.7 707.7 712.7 717.7 722.7 752.7 747.7 742.7 737.7 732.7 727.7 -o k19_vertical_scan_cluster_position
