#!/bin/bash
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



files_high=(/scratch/data/testbeam201905/Data_runs/Run_20190505_100835.dat)

calibration=/scratch/data/testbeam201905/Data_runs/Calibration_20190504_230648.dat.ymlcalib.root
calibration_high=/scratch/data/testbeam201905/Data_runs/Calibration_20190505_080154.dat.ymlcalib.root

pedestal_ext=".tree_pedestal.root"


for i in ${files[*]}
do
	./bin/pedestal_tree $i $calibration
	#./bin/analysisExternal_tree $i $calibration $i$pedestal_ext
	
done

#for i in ${files_high[*]}
#do
	#./bin/pedestal_tree $i $calibration_high
	#./bin/analysisExternal_tree $i $calibration_high $i$pedestal_ext
	
#done

#source ./scripts/tb201905/ext_plots_script.sh
