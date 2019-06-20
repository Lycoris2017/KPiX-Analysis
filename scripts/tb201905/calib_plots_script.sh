#!/bin/bash
#calib_files_normal_gain=(/scratch/data/testbeam201905/Data_runs/Calibration_20190504_230648.dat /scratch/data/testbeam201905/debug_runs/Calibration_20190513_152812.dat)

##/scratch/data/testbeam201905/Calibration_20190502_083214.dat /scratch/data/testbeam201905/Calibration_20190502_180620.dat /scratch/data/testbeam201905/Calibration_20190503_121847.dat /scratch/data/testbeam201905/Calibration_20190504_111600.dat


ymlcalib_extension=.ymlcalib.root
##

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
#for i in ${kpix[*]}
#do
	##python2.7 python/plot_producer.py /scratch/data/testbeam201905/Data_runs/Calibration_20190523_095745.dat.ymlcalib.root /scratch/data/testbeam201905/Data_runs/Calibration_20190523_115002.dat.ymlcalib.root -n "slope_0RMS_"$i"_b0" -d 'hist e same' -o "calibration_swapping_RMS_test_"$i --legend Normal Swapped
	#python2.7 python/plot_producer.py /scratch/data/testbeam201905/Data_runs/Calibration_20190504_230648.dat.ymlcalib.root -n "slope_"$i -d "h e same" -o "bucket_compare_normalgain_"$i --legend b0 b1 --refuse b2 b3

#done

## To show dieter that the CalDacCount has an impact
calib_debug_files=(/scratch/data/testbeam201905/debug_runs/Calibration_20190613_155332.dat /scratch/data/testbeam201905/debug_runs/Calibration_20190613_165816.dat)
same=""
for i in ${calib_debug_files[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	same="$same $i$ymlcalib_extension"
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestalsRMS_fc_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare_downstream -b 0 --refuse 10 11 --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2 --xrange -0.5 4
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestalsRMS_fc_k -d "h e same" -k 6 7 8 9 10 11 -o pedestalRMS_compare_upstream -b 0 --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2 --xrange -0.5 4
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_k -d "h e same" -k 0 1 2 3 4 5 -o slope_compare_downstream -b 0 --refuse 10 11 --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_k -d "h e same" -k 6 7 8 9 10 11 -o slope_compare_upstream -b 0 --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2
	
	#python python/plot_producer.py $i$ymlcalib_extension -n pearson_correlation -d "h e same" -k 0 1 2 3 4 5 -o pearson_compare_downstream --xrange -1 2 -b 0 --refuse 10 11 --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2
	#python python/plot_producer.py $i$ymlcalib_extension -n pearson_correlation -d "h e same" -k 6 7 8 9 10 11 -o pearson_compare_upstream --xrange -1 2 -b 0 --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2

	python python/plot_producer.py $i$ymlcalib_extension -n calib_fluctuate_2D -d "colz" --yrange 200 255 --xrange 0 8191 --zlog True

done
#for i in ${kpix[*]}
#do
	#python python/plot_producer.py $same -n "slope_"$i"_" -d "h e same" -b 0 -o "CalDacCount_compare_"$i --legend CalDacCount1 CalDacCount5
	#python python/plot_producer.py $same -n "slope_"$i"_" -d "h e same" -b 0 -o "CalDacCount_compare_-5to5_"$i --legend CalDacCount1 CalDacCount5 --xrange -5 5
	#python python/plot_producer.py $same -n "slope_"$i"_" -d "h e same" -b 0 -o "CalDacCount_compare_5to35_"$i --legend CalDacCount1 CalDacCount5 --xrange 5 35
#done

##To look into the lower cntrlHoldTime settings in calibration
#calib_cntrlHoldTime=(/scratch/data/testbeam201905/debug_runs/Calibration_20190611_134732.dat /scratch/data/testbeam201905/debug_runs/Calibration_20190611_152257.dat /scratch/data/testbeam201905/debug_runs/Calibration_20190611_170251.dat /scratch/data/testbeam201905/debug_runs/Calibration_20190613_165816.dat /scratch/data/testbeam201905/debug_runs/Calibration_20190614_125706.dat)
#same=""
#for i in ${calib_cntrlHoldTime[*]}
#do
	#substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	#same="$same $i$ymlcalib_extension"
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestalsRMS_fc_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare_downstream -b 0 --refuse 10 11 --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2 --xrange -0.5 4
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestalsRMS_fc_k -d "h e same" -k 6 7 8 9 10 11 -o pedestalRMS_compare_upstream -b 0 --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2 --xrange -0.5 4
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_k -d "h e same" -k 0 1 2 3 4 5 -o slope_compare_downstream -b 0 --refuse 10 11 --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_k -d "h e same" -k 6 7 8 9 10 11 -o slope_compare_upstream -b 0 --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2
	
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestals_fc_k -d "h e same" -k 6 7 8 9 10 11 -o pedestalsfC_compare_upstream -b 0 --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestals_k -d "h e same" -k 6 7 8 9 10 11 -o pedestals_compare_upstream -b 0 --legend S34-1 S34-2 S53-1 S53-2 S55-1 S55-2
	
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestals_fc_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare_downstream -b 0 --refuse 10 11 --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestals_k -d "h e same" -k 0 1 2 3 4 5 -o pedestalRMS_compare_downstream -b 0 --refuse 10 11 --legend S59-1 S59-2 S43-1 S43-2 S1st-1 S1st-2
	
#done
#for i in ${kpix[*]}
#do
	#python python/plot_producer.py $same -n "slope_"$i"_" -d "h e same" -b 0 -o "CntrHoldTime_compare_slope_"$i --legend 64x_320ns 24x_320ns 24x_160ns 64x_320ns 32x_320ns
	#python python/plot_producer.py $same -n "pedestals_fc_"$i -d "h e same"  -o "CntrHoldTime_compare_pedestalsfC_"$i -b 0 --legend 64x_320ns 24x_320ns 24x_160ns 64x_320ns 32x_320ns
	#python python/plot_producer.py $same -n "pedestals_"$i -d "h e same"  -o "CntrHoldTime_compare_pedestals_0to2000_"$i -b 0 --legend 64x_320ns 24x_320ns 24x_160ns 64x_320ns 32x_320ns --xrange 0 2000
	#python python/plot_producer.py $same -n "pedestals_"$i -d "h e same"  -o "CntrHoldTime_compare_pedestals_"$i -b 0 --legend 64x_320ns 24x_320ns 24x_160ns 64x_320ns 32x_320ns
#done
