#!/bin/bash
calib_files20190402_normal_gain=(/scratch/data/2019_04_02_16_45_04.bin /scratch/data/2019_04_02_17_09_58.bin)

calib_files20190403_normal_gain=(/scratch/data/2019_04_03_14_20_27.bin /scratch/data/2019_04_03_14_44_39.bin)
calib_files20190403_high_gain=(/scratch/data/2019_04_03_15_11_44.bin /scratch/data/2019_04_03_15_42_04.bin)


calib_files20190404_normal_gain=(/scratch/data/2019_04_04_14_37_03.bin /scratch/data/2019_04_04_15_13_41.bin) 
calib_files20190404_normal_gain512=(/scratch/data/2019_04_04_15_01_06.bin /scratch/data/2019_04_04_15_37_00.bin)

calib_newDAQ_20190404_05=(/home/lycoris-dev/workspace/kpix/software/data_20190404_155926.dat /home/lycoris-dev/workspace/kpix/software/data_20190404_161542.dat /home/lycoris-dev/workspace/kpix/software/data_20190404_162425.dat /scratch/data/newdaq/data_20190405_153930.dat)
calib_newDAQ_20190408=(/scratch/data/newdaq/data_20190408_133710.dat /scratch/data/newdaq/data_20190408_135420.dat /scratch/data/newdaq/data_20190408_150832.dat)
calib_files20190408_normal_gain=(/scratch/data/2019_04_08_13_05_13.bin)

newdaq_olddaq=(/scratch/data/newdaq/data_20190408_150832.dat /scratch/data/2019_04_08_13_05_13.bin)


calib_extension=.newCalib.root
ymlcalib_extension=.ymlcalib.root
#same=""

#for i in ${calib_files20190402_normal_gain[*]}
#do
	#substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	#same="$same $i$calib_extension"
#done

#python python/plot_producer.py $same -n pedestals_fc_k -b 0 -d 'hist e same' --legend S41-2 S41-1 -o S41_pedestal_comparison --rebin 8
#python python/plot_producer.py $same -n slope_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_comparison
#python python/plot_producer.py $same -n slope_0_127_k -b 0  -d 'hist e same'   --legend S41-2 S41-1 -o S41_slope_0_127_comparison
#python python/plot_producer.py $same -n RMS_fc_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_RMSfc_comparison --rebin 8	
#python python/plot_producer.py $same -n RMS_fc_0_127_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_RMSfc_0_127_comparison --rebin 8
#python python/plot_producer.py $same -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_vs_channel_comparison	
#python python/plot_producer.py $same -n slope_vs_right_strip_k  -b 0 -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_vs_right_strip_comparison	

#same=""

#for i in ${calib_files20190403_normal_gain[*]}
#do
	#substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	#same="$same $i$calib_extension"
#done

#python python/plot_producer.py $same -n pedestals_fc_k -b 0 -d 'hist e same' --legend S41-2 S41-1 -o S41_pedestal_comparison --rebin 8
#python python/plot_producer.py $same -n slope_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_comparison
#python python/plot_producer.py $same -n RMS_fc_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_RMSfc_comparison --rebin 8	
#python python/plot_producer.py $same -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_vs_channel_comparison	

#same=""

#for i in ${calib_files20190403_high_gain[*]}
#do
	#substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	#same="$same $i$calib_extension"
#done

#python python/plot_producer.py $same -n pedestals_fc_k -b 0 -d 'hist e same' --legend S41-2 S41-1 -o S41_pedestal_comparison --rebin 8
#python python/plot_producer.py $same -n slope_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_comparison
#python python/plot_producer.py $same -n RMS_fc_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_RMSfc_comparison --rebin 8	
#python python/plot_producer.py $same -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_vs_channel_comparison	

#same=""

#for i in ${calib_files20190404_normal_gain[*]}
#do
	#substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	#same="$same $i$calib_extension"
#done

#python python/plot_producer.py $same -n pedestals_fc_k -b 0 -d 'hist e same' --legend S41-2 S41-1 -o S41_pedestal_comparison --rebin 8
#python python/plot_producer.py $same -n slope_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_comparison
#python python/plot_producer.py $same -n RMS_fc_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_RMSfc_comparison --rebin 8	
#python python/plot_producer.py $same -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_vs_channel_comparison	

#same=""

#for i in ${calib_files20190404_normal_gain512[*]}
#do
	#substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	#same="$same $i$calib_extension"
#done

#python python/plot_producer.py $same -n pedestals_fc_k -b 0 -d 'hist e same' --legend S41-2 S41-1 -o S41_pedestal_comparison --rebin 8
#python python/plot_producer.py $same -n slope_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_comparison
#python python/plot_producer.py $same -n RMS_fc_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_RMSfc_comparison --rebin 8	
#python python/plot_producer.py $same -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_vs_channel_comparison	

#same=""

#for i in ${calib_newDAQ_20190404_05[*]}
#do
	
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestals_fc_k -b 0 -d 'hist e same' --legend S41-1 S41-2 -o S41_pedestal_comparison --rebin 8
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_k -b 0  -d 'hist e same'  --legend S41-1 S41-2 -o S41_slope_comparison
	#python python/plot_producer.py $i$ymlcalib_extension -n RMS_fc_k -b 0  -d 'hist e same'  --legend S41-1 S41-2 -o S41_RMSfc_comparison --rebin 8	
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend S41-1 S41-2 -o S41_slope_vs_channel_comparison	
	#python python/plot_producer.py $i$ymlcalib_extension -n calib_ -c 0117 0390 0417 0512 0700 0892  -b 0 -d '' 
#done

for i in ${calib_newDAQ_20190408[*]}
do
	
	#python python/plot_producer.py $i$ymlcalib_extension -n pedestals_fc_k -b 0 -d 'hist e same' --legend S41-1 S41-2 -o S41_pedestal_comparison --rebin 8
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_k -b 0  -d 'hist e same'  --legend S41-1 S41-2 -o S41_slope_comparison
	#python python/plot_producer.py $i$ymlcalib_extension -n RMS_fc_k -b 0  -d 'hist e same'  --legend S41-1 S41-2 -o S41_RMSfc_comparison --rebin 8	
	#python python/plot_producer.py $i$ymlcalib_extension -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend S41-1 S41-2 -o S41_slope_vs_channel_comparison	
	python python/plot_producer.py $i$calib_extension -n calib_ -c 0001 0053 0078 0117 0390 0417 0512 0700 0892  -b 0 -d ''
done

for i in ${calib_files20190408_normal_gain[*]}
do
	#python python/plot_producer.py $i$calib_extension -n pedestals_fc_k -b 0 -d 'hist e same' --legend S41-2 -o S41-2_pedestal_comparison --rebin 8
	#python python/plot_producer.py $i$calib_extension -n slope_k -b 0  -d 'hist e same'  --legend S41-2 -o S41-2_slope_comparison
	#python python/plot_producer.py $i$calib_extension -n RMS_fc_k -b 0  -d 'hist e same'  --legend S41-2 -o S41-2_RMSfc_comparison --rebin 8	
	#python python/plot_producer.py $i$calib_extension -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend S41-2 -o S41-2_slope_vs_channel_comparison	
	python python/plot_producer.py $i$calib_extension -n calib_ -c 0001 0053 0078 0117 0390 0417 0512 0700 0892  -b 0 -d '' 
done

#python python/plot_producer.py $same -n pedestals_fc_k -b 0 -d 'hist e same' --legend S41-2 S41-1 -o S41_pedestal_comparison --rebin 8
#python python/plot_producer.py $same -n slope_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_comparison
#python python/plot_producer.py $same -n RMS_fc_k -b 0  -d 'hist e same'  --legend S41-2 S41-1 -o S41_RMSfc_comparison --rebin 8	
#python python/plot_producer.py $same -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend S41-2 S41-1 -o S41_slope_vs_channel_comparison	


#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n pedestals_fc_k -b 0 -d 'hist e same' --legend newDAQ oldDAQ -o DAQ_comparison_pedestal_comparison --rebin 8
#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n slope_k -b 0  -d 'hist e same'  --legend newDAQ oldDAQ -o DAQ_comparison_slope_comparison
#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n RMS_fc_k -b 0  -d 'hist e same'  --legend newDAQ oldDAQ -o DAQ_comparison_RMSfc_comparison --rebin 8	
#python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newCalib.root -n slope_vs_channel_k  -b 0 -d 'hist e same'  --legend newDAQ oldDAQ -o DAQ_comparison_slope_vs_channel_comparison	
##python python/plot_producer.py /scratch/data/newdaq/data_20190408_150832.dat.ymlcalib.root /scratch/data/2019_04_08_13_05_13.bin.newcalib.root -n calib_ -c 0117 0390 0417 0512 0700 0892  -b 0 -d ''
