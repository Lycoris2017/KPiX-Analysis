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
/scratch/data/testbeam201905/Data_runs/Run_20190505_172140.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_175332.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_181637.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_182346.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_183918.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_183918.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_190055.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_192134.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_194638.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_202335.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_195639.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_204113.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_210404.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_211206.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_213817.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_215704.dat
/scratch/data/testbeam201905/Data_runs/Run_20190505_220458.dat)

calibration=/scratch/data/testbeam201905/Data_runs/Calibration_20190504_230648.dat.ymlcalib.root



for i in ${files[*]}
do
	./bin/analysis_newdaq $i
	
done

