#!/bin/bash
large_file_70v=/scratch/data/testbeam201902/2019_02_13_15_37_39.bin
bias_scan_files1=(/scratch/data/testbeam201902/2019_02_13_13_46_28.bin /scratch/data/testbeam201902/2019_02_13_13_50_24.bin /scratch/data/testbeam201902/2019_02_13_13_57_50.bin /scratch/data/testbeam201902/2019_02_13_14_01_29.bin /scratch/data/testbeam201902/2019_02_13_14_09_01.bin) 
bias_scan_files2=(/scratch/data/testbeam201902/2019_02_13_14_13_53.bin /scratch/data/testbeam201902/2019_02_13_14_26_12.bin /scratch/data/testbeam201902/2019_02_13_14_30_39.bin /scratch/data/testbeam201902/2019_02_13_14_40_41.bin /scratch/data/testbeam201902/2019_02_13_14_45_36.bin) 
bias_scan_files3=(/scratch/data/testbeam201902/2019_02_13_14_57_20.bin /scratch/data/testbeam201902/2019_02_13_15_15_52.bin /scratch/data/testbeam201902/2019_02_13_15_26_42.bin /scratch/data/testbeam201902/2019_02_13_15_31_42.bin) 

bias_scan_files_sample=(/scratch/data/testbeam201902/2019_02_13_13_50_24.bin /scratch/data/testbeam201902/2019_02_13_14_13_53.bin /scratch/data/testbeam201902/2019_02_13_14_45_36.bin /scratch/data/testbeam201902/2019_02_13_15_31_42.bin)



extension=.pedestal.external.root
same=""
same1=""
same2=""
same3=""
same_sample=""


for i in ${bias_scan_files1[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	same="$same $i"_"$substring$extension"
	same1="$same1 $i"_"$substring$extension"
	
done

for i in ${bias_scan_files2[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	same="$same $i"_"$substring$extension"
	same2="$same2 $i"_"$substring$extension"
done

for i in ${bias_scan_files3[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	same="$same $i"_"$substring$extension"
	same3="$same3 $i"_"$substring$extension"
done

for i in ${bias_scan_files_sample[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	same_sample="$same_sample $i"_"$substring$extension"
done



#python2.7 python/plot_producer.py $same -n cluster_size k17 -b0  -d 'hist e same' --legend 70V 80V 90V 100V 110V 120V 130V 140V 150V 160V 170V 180V 190V 200V -o cluster_size_k17_bias_scan_70to200
#python2.7 python/plot_producer.py $same -n cluster_size k19 -b0  -d 'hist e same' --legend 70V 80V 90V 100V 110V 120V 130V 140V 150V 160V 170V 180V 190V 200V -o cluster_size_k19_bias_scan_70to200


#python2.7 python/plot_producer.py $same1 -n cluster_size k17 -b0  -d 'hist e same' --legend 70V 80V 90V 100V 110V  -o cluster_size_k17_bias_scan_70to110
#python2.7 python/plot_producer.py $same1 -n cluster_size k19 -b0  -d 'hist e same' --legend 70V 80V 90V 100V 110V  -o cluster_size_k19_bias_scan_70to110


#python2.7 python/plot_producer.py $same2 -n cluster_size k17 -b0  -d 'hist e same' --legend 120V 130V 140V 150V 160V  -o cluster_size_k17_bias_scan_120to160
#python2.7 python/plot_producer.py $same2 -n cluster_size k19 -b0  -d 'hist e same' --legend 120V 130V 140V 150V 160V  -o cluster_size_k19_bias_scan_120to160


#python2.7 python/plot_producer.py $same3 -n cluster_size k17 -b0  -d 'hist e same' --legend 170V 180V 190V 200V -o cluster_size_k17_bias_scan_170to200
#python2.7 python/plot_producer.py $same3 -n cluster_size k19 -b0  -d 'hist e same' --legend 170V 180V 190V 200V -o cluster_size_k19_bias_scan_170to200


python2.7 python/plot_producer.py $same_sample -n cluster_size k17 -b0  -d 'hist e same' --legend 80V 120V 160V 200V -o cluster_size_k17_bias_scan_sample
python2.7 python/plot_producer.py $same_sample -n cluster_size k19 -b0  -d 'hist e same' --legend 80V 120V 160V 200V -o cluster_size_k19_bias_scan_sample

python2.7 python/plot_producer.py $same_sample -n median_made k17 -b0  -d 'hist e same' --xrange -3 3 --legend 80V 120V 160V 200V -o cluster_size_k17_bias_scan_noise
python2.7 python/plot_producer.py $same_sample -n median_made k19 -b0  -d 'hist e same' --xrange -3 3 --legend 80V 120V 160V 200V -o cluster_size_k19_bias_scan_noise
