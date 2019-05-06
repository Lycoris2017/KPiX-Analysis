#!/bin/bash
files=(/scratch/data/testbeam201905/Data_runs/Run_20190504_225322.dat /scratch/data/testbeam201905/Data_runs/Run_20190505_094624.dat /scratch/data/testbeam201905/Data_runs/Run_20190504_225553.dat /scratch/data/testbeam201905/Data_runs/Run_20190504_225931.dat)
files_high=(/scratch/data/testbeam201905/Data_runs/Run_20190505_100835.dat)

calibration=/scratch/data/testbeam201905/Data_runs/Calibration_20190504_230648.dat.ymlcalib.root
calibration_high=/scratch/data/testbeam201905/Data_runs/Calibration_20190505_080154.dat.ymlcalib.root

pedestal_ext=".pedestal.root"


for i in ${files[*]}
do
	./bin/pedestal $i $calibration
	./bin/analysisExternal $i $calibration $i$pedestal_ext
	
done

for i in ${files_high[*]}
do
	./bin/pedestal $i $calibration_high
	./bin/analysisExternal $i $calibration_high $i$pedestal_ext
	
done
