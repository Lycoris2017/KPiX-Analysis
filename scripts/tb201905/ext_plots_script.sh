#!/bin/bash
#files0215=(/scratch/data/testbeam201902/2019_02_15_11_22_33.bin /scratch/data/testbeam201902/2019_02_15_11_26_40.bin /scratch/data/testbeam201902/2019_02_15_11_30_22.bin /scratch/data/testbeam201902/2019_02_15_11_32_50.bin)
#files0213=(/scratch/data/testbeam201902/2019_02_13_13_42_33.bin /scratch/data/testbeam201902/2019_02_13_13_46_28.bin  /scratch/data/testbeam201902/2019_02_13_15_37_39.bin /scratch/data/testbeam201902/2019_02_13_13_39_57.bin /scratch/data/testbeam201902/2019_02_13_13_36_34.bin /scratch/data/testbeam201902/2019_02_13_13_33_12.bin /scratch/data/testbeam201902/2019_02_13_13_28_11.bin /scratch/data/testbeam201902/2019_02_13_13_24_39.bin /scratch/data/testbeam201902/2019_02_13_13_19_57.bin /scratch/data/testbeam201902/2019_02_13_12_14_12.bin /scratch/data/testbeam201902/2019_02_13_12_10_35.bin /scratch/data/testbeam201902/2019_02_13_12_00_47.bin /scratch/data/testbeam201902/2019_02_13_11_56_38.bin)

files=(/scratch/data/testbeam201905/Data_runs/Run_20190504_225322.dat
/scratch/data/testbeam201905/Data_runs/Run_20190504_225553.dat
/scratch/data/testbeam201905/Data_runs/Run_20190504_225931.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_094624.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_111442.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_113435.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_121902.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_122701.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_125110.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_130524.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_133249.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_135207.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_135952.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_143637.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_145353.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_150111.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_153442.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_154212.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_160444.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_161533.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_163916.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_164842.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_171322.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_172837.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_180052.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_180652.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_183051.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_183051.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_184618.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_185353.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_192841.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_193859.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_201615.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_200442.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_204838.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_205604.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_211922.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_212658.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_214958.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_221221.dat
)


file2=/scratch/data/testbeam201905/Data_runs/Run_20190504_225553.dat
file3=/scratch/data/testbeam201905/Data_runs/Run_20190504_225931.dat

#pedestal0215=2019_02_15_11_22_33


extension=.pedestal.external.root

substring=$(echo $file2| cut -d'/' -f 6 | cut -d '.' -f 1)
python2.7 python/plot_producer.py $file2"_"$substring$extension -n fc_response_median_made _k -k 0 3 5 7 8 10 -b 0 -d 'hist e same' -o fc_response_median_made_compare_upperRow --refuse 11 9 16 12 13 _c --legend S59-1 S43-2 S1st-2 S34-2 S53-1 S55-1
substring=$(echo $file3| cut -d'/' -f 6 | cut -d '.' -f 1)
python2.7 python/plot_producer.py $file3"_"$substring$extension -n fc_response_median_made _k -k 1 2 4 6 9 11  -b 0 -d 'hist e same' -o fc_response_median_made_compare_lowerRow --refuse 10 8 7 16 13 12 _c --legend S59-2 S43-1 S1st-1 S34-1 S53-2 S55-2

for i in ${files[*]}
do
	substring=$(echo $i| cut -d'/' -f 6 | cut -d '.' -f 1)
	#echo $substring
	#echo $i"_"$substring$extension
	#python2.7 python/plot_producer.py $i"_"$substring$extension -n slopes_vs_stripr -b 0 -d 'hist e'
	#python2.7 python/plot_producer.py $i"_"$substring$extension -n slopes_vs_channel -b 0 -d 'hist e'

	#python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n fc_response_median_made _k -b 0 -d 'hist e' --refuse 28 30
	
	
	python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_median_made _k  -k 0 1 2 3 4 5 -b 0 -d 'hist e same' -o fc_response_median_made_compare_downstream --refuse 10 11 16 12 13 _c --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2 
	python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_median_made _k  -k 6 7 8 9 10 11 -b 0 -d 'hist e same' -o fc_response_median_made_compare_upstream --refuse 16 13 12 _c --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2
	
	python2.7 python/plot_producer.py $i"_"$substring$extension -n noise_distribution _b0  -k 0 1 2 3 4 5 -b 0 -d 'hist e same' -o noise_distribution_downstream --refuse 10 11 16 12 13 _c --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2 
	python2.7 python/plot_producer.py $i"_"$substring$extension -n noise_distribution _b0  -k 6 7 8 9 10 11 -b 0 -d 'hist e same' -o noise_distribution_compare_upstream --refuse 16 13 12 _c --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2
	#python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_median_made _k  -k 0 3 5 7 8 10 -b 0 -d 'hist e same' -o fc_response_median_made_compare_upperRow --refuse  _c --legend S59-1 S43-2 S1st-2 S34-2 S53-1 S55-1
	#python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_median_made _k  -k 1 2 4 6 9 11 -b 0 -d 'hist e same' -o fc_response_median_made_compare_lowerRow --refuse  _c --legend S59-2 S43-1 S1st-1 S34-1 S53-2 S55-2
	#python2.7 python/plot_producer.py $i"_"$substring$extension -n fc_response_median_made _k -b 0 -d 'hist e'
	
	
	#python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n cluster_position -b 0 -d 'h e' --refuse 28 30
	#python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n cluster offset -d 'h e' --refuse 28 30
	
	#python2.7 python/plot_producer.py $i"_"$pedestal0215$extension -n cluster_charge -b 0 -d 'h e' --refuse 28 30
	
done
