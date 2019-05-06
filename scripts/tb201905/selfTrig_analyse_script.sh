#!/bin/bash
files=(/scratch/data/testbeam201905/Data_runs/Run_20190505_115706.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_121155.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_123403.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_124320.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_131633.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_132512.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_134141.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_140910.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_142913.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_144632.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_150822.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_151603.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_152740.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_154909.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_155731.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_162237.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_163133.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_165605.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_170528.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_172140.dat)

calibration=/scratch/data/testbeam201905/Data_runs/Calibration_20190504_230648.dat.ymlcalib.root



for i in ${files[*]}
do
	./bin/analysis_newdaq $i
	
done

