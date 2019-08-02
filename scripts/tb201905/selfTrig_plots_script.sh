#!/bin/bash
files=(/scratch/data/testbeam201905/Data_runs/Run_20190505_142913.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_140910.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_134141.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_132512.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_131633.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_124320.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_123403.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_121155.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_115706.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_165605.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_144632.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_150822.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_151603.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_154909.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_152740.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_155731.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_162237.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_163133.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_170528.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_172140.dat)


#


extension=.root

output_name=(strip_entries_447_221 
strip_entries_437_221 
strip_entries_427_221 
strip_entries_417_221 
strip_entries_407_221 
strip_entries_397_221 
strip_entries_387_221 
strip_entries_382_221 
strip_entries_377_221 
strip_entries_377_221 
strip_entries_372_221 
strip_entries_367_221 
strip_entries_362_221 
strip_entries_362_221 
strip_entries_352_221 
strip_entries_342_221 
strip_entries_332_221 
strip_entries_322_221 
strip_entries_312_221
strip_entries_302_221)

count=0

for i in ${files[*]}
do
	
	python2.7 python/plot_producer.py $i$extension -n strip_entries _k  -k 0 3 5 7 8 10 -b 0 -d 'hist e same' -o ${output_name[$count]} --refuse timed --legend S59-1 S43-2 S1st-2 S34-2 S53-1 S55-1
	python2.7 python/plot_producer.py $i$extension -n timed_ strip_entries _k  -k 0 3 5 7 8 10 -b 0 -d 'hist e same' -o timed_${output_name[$count]} --legend S59-1 S43-2 S1st-2 S34-2 S53-1 S55-1
	count=$[count + 1]
done
