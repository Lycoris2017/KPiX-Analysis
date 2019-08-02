#!bin/sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib

## External Trigger w/ TS
# /scratch/data/tracker_test/eudaq//testbeam_201808/2018_08_18_20_06_16.bin
## Auto Trigger w/ TS
# /scratch/data/tracker_test/2018_08_19_22_32_03.bin
## Auto Trigger


#bin/count /scratch/data/tracker_test/eudaq//testbeam_201808/2018_08_18_20_06_16.bin
#./bin/ntupleTest /scratch/data/tracker_test/eudaq//testbeam_201808/2018_08_18_20_06_16.bin
./bin/ntupleTest /scratch/data/tracker_test/2018_08_19_22_32_03.bin
