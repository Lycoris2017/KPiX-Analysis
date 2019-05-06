#!/bin/bash
calib_files20190402_normal_gain=(/scratch/data/testbeam201905/Calibration_20190504_230648.dat)
calib_files20190402_normal_gain=(/scratch/data/testbeam201905/Calibration_20190504_230648.dat)
#/scratch/data/testbeam201905/Calibration_20190502_083214.dat /scratch/data/testbeam201905/Calibration_20190502_180620.dat /scratch/data/testbeam201905/Calibration_20190503_121847.dat /scratch/data/testbeam201905/Calibration_20190504_111600.dat


ymlcalib_extension=.ymlcalib.root
#same=""

for i in ${calib_files20190402_normal_gain[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	same="$same $i$ymlcalib_extension"
	python python/plot_producer.py $i$ymlcalib_extension -n pedestalsRMS_fc_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare -b 0 --refuse 10 11 --legend k0 k1 k2 k3 k4 k5 --xrange -0.5 4
	python python/plot_producer.py $i$ymlcalib_extension -n pedestalsRMS_fc_k -d "h e same" -k 6 7 8 9 10 11 -o pedestalRMS_compare2 -b 0 --legend k6 k7 k8 k9 k10 k11 --xrange -0.5 4
	python python/plot_producer.py $i$ymlcalib_extension -n slope_k -d "h e same" -k 0 1 2 3 4 5 -o slope_compare -b 0 --refuse 10 11 --legend k0 k1 k2 k3 k4 k5
	python python/plot_producer.py $i$ymlcalib_extension -n slope_k -d "h e same" -k 6 7 8 9 10 11 -o slope_compare2 -b 0 --legend k6 k7 k8 k9 k10 k11
	
done

#python python/plot_producer.py $same -n pedestalsRMS_fc_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare -b 0 --refuse 10 11 --legend k0 k1 k2 k3 k4 k5 --xrange -0.5 4
#python python/plot_producer.py $same -n pedestalsRMS_fc_k -d "h e same" -k 6 7 8 9 10 11 -o pedestalRMS_compare2 -b 0 --legend k6 k7 k8 k9 k10 k11 --xrange -0.5 4
#python python/plot_producer.py $same -n slope_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare -b 0 --refuse 10 11 --legend k0 k1 k2 k3 k4 k5 --xrange -0.5 4
#python python/plot_producer.py $same -n slope_k -d "h e same" -k 6 7 8 9 10 11 -o pedestalRMS_compare2 -b 0 --legend k6 k7 k8 k9 k10 k11 --xrange -0.5 4








#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n pedestals_fc_k -b 0 -d 'hist e same' --legend newDAQ oldDAQ -o DAQ_comparison_pedestal_comparison --rebin 8
#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n slope_k -b 0  -d 'hist e same'  --legend newDAQ oldDAQ -o DAQ_comparison_slope_comparison
#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n RMS_fc_k -b 0  -d 'hist e same'  --legend newDAQ oldDAQ -o DAQ_comparison_RMSfc_comparison --rebin 8	
#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend newDAQ oldDAQ -o DAQ_comparison_slope_vs_channel_comparison	
##python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newcalib.root -n calib_ -c 0117 0390 0417 0512 0700 0892  -b 0 -d ''
