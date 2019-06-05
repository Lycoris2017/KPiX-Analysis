#!/bin/bash
#calib_files_normal_gain=(/scratch/data/testbeam201905/Data_runs/Calibration_20190504_230648.dat /scratch/data/testbeam201905/debug_runs/Calibration_20190513_152812.dat)

##/scratch/data/testbeam201905/Calibration_20190502_083214.dat /scratch/data/testbeam201905/Calibration_20190502_180620.dat /scratch/data/testbeam201905/Calibration_20190503_121847.dat /scratch/data/testbeam201905/Calibration_20190504_111600.dat


#ymlcalib_extension=.ymlcalib.root
##same=""

#for i in ${calib_files_normal_gain[*]}
#do
	#substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	#same="$same $i$ymlcalib_extension"
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestalsRMS_fc_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare_downstream -b 0 --refuse 10 11 --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2 --xrange -0.5 4
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestalsRMS_fc_k -d "h e same" -k 6 7 8 9 10 11 -o pedestalRMS_compare_upstream -b 0 --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2 --xrange -0.5 4
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_k -d "h e same" -k 0 1 2 3 4 5 -o slope_compare_downstream -b 0 --refuse 10 11 --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_k -d "h e same" -k 6 7 8 9 10 11 -o slope_compare_upstream -b 0 --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2
	
#done

#python python/plot_producer.py $same -n pedestalsRMS_fc_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare -b 0 --refuse 10 11 --legend k0 k1 k2 k3 k4 k5 --xrange -0.5 4
#python python/plot_producer.py $same -n pedestalsRMS_fc_k -d "h e same" -k 6 7 8 9 10 11 -o pedestalRMS_compare2 -b 0 --legend k6 k7 k8 k9 k10 k11 --xrange -0.5 4
#python python/plot_producer.py $same -n slope_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare -b 0 --refuse 10 11 --legend k0 k1 k2 k3 k4 k5 --xrange -0.5 4
#python python/plot_producer.py $same -n slope_k -d "h e same" -k 6 7 8 9 10 11 -o pedestalRMS_compare2 -b 0 --legend k6 k7 k8 k9 k10 k11 --xrange -0.5 4








#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n pedestals_fc_k -b 0 -d 'hist e same' --legend newDAQ oldDAQ -o DAQ_comparison_pedestal_comparison --rebin 8
#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n slope_k -b 0  -d 'hist e same'  --legend newDAQ oldDAQ -o DAQ_comparison_slope_comparison
#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n RMS_fc_k -b 0  -d 'hist e same'  --legend newDAQ oldDAQ -o DAQ_comparison_RMSfc_comparison --rebin 8	
#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend newDAQ oldDAQ -o DAQ_comparison_slope_vs_channel_comparison	
##python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newcalib.root -n calib_ -c 0117 0390 0417 0512 0700 0892  -b 0 -d ''


kpix=(k0 k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11)
for i in ${kpix[*]}
do
	python2.7 python/plot_producer.py /scratch/data/testbeam201905/Data_runs/Calibration_20190523_095745.dat.ymlcalib.root /scratch/data/testbeam201905/Data_runs/Calibration_20190523_115002.dat.ymlcalib.root -n "slope_0RMS_"$i"_b0" -d 'hist e same' -o "calibration_swapping_RMS_test_"$i --legend Normal Swapped
done
