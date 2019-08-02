#!/bin/bash

calib_ext=".newCalib.root"
pedestal_ext=".pedestal.root"

data_file=$1
calib_file=$2
xml_file=$3


#echo "$data_file"

./bin/new_CalibrationFitter $xml_file $calib_file
./bin/pedestal $data_file $calib_file$calib_ext
./bin/analysisExternal $data_file $calib_file$calib_ext $data_file$pedestal_ext

