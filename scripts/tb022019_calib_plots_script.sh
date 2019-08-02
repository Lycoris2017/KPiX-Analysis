#!/bin/bash
calib_files=(/scratch/data/testbeam201902/2019_02_11_10_43_00.bin /scratch/data/testbeam201902/2019_02_12_17_19_18.bin /scratch/data/testbeam201902/2019_02_13_12_21_15.bin /scratch/data/testbeam201902/2019_02_14_11_54_29.bin /scratch/data/testbeam201902/2019_02_15_10_03_05.bin /scratch/data/testbeam201902/2019_02_20_10_41_07.bin)
calib_extension=.newCalib.root
kpix=(17 19 22 28 30)

same=""
dates=""

for i in ${calib_files[*]}
do
	substring=$(echo $i| cut -d'/' -f 5 | cut -d '.' -f 1)
	same="$same $i$calib_extension"
	dates="$dates $substring"
	##for k in ${kpix[*]}
	##do
	#python python/plot_producer.py $i$calib_extension -n pedestals_fc_k -b 0 -d 'hist e' --rebin 8
	#python python/plot_producer.py $i$calib_extension -n slope_k -b 0  -d 'hist e'
	#python python/plot_producer.py $i$calib_extension -n slope_0_127_k -b 0  -d 'hist e'
	#python python/plot_producer.py $i$calib_extension -n RMS_fc_k  -b 0 -d 'hist e' --rebin 8
	#python python/plot_producer.py $i$calib_extension -n RMS_fc_0_127_k  -b 0 -d 'hist e' --rebin 8
	#python python/plot_producer.py $i$calib_extension -n RMS_fc _disc_ -b 0 -d 'hist e' --rebin 8
	#python python/plot_producer.py $i$calib_extension -n RMS_fc _conn_ -b 0 -d 'hist e' --rebin 8
	#python python/plot_producer.py $i$calib_extension -n slope_vs_channel_k  -b 0 -d 'hist e'
	#python python/plot_producer.py $i$calib_extension -n slope_vs_right_strip_k  -b 0 -d 'hist e'
	
	#python python/plot_producer.py $i$calib_extension -n hist_10 -b 0 -c 0211 0332 0467 0784 0087 1000 -d 'hist e'
	
	#python python/plot_producer.py $i$calib_extension -n calib_10 -b 0 -c 0211 0332 0467 0784 0087 1000 -d 'e'
	
	##done	
	#python python/plot_producer.py $i$calib_extension -n pedestals_fc_k -b 0 -d 'hist e same' --refuse 28 30 --legend S59-2 S43-2 -o tb201902_pedestal_comparison.png --rebin 8
	#python python/plot_producer.py $i$calib_extension -n slope_k -b 0  -d 'hist e same' --refuse 28 30  --legend S59-2 S43-2 -o tb201902_slope_comparison.png
	#python python/plot_producer.py $i$calib_extension -n slope_0_127_k -b 0  -d 'hist e same' --refuse 28 30  --legend S59-2 S43-2 -o tb201902_slope_0_127_comparison.png
	#python python/plot_producer.py $i$calib_extension -n RMS_fc_k -b 0  -d 'hist e same' --refuse 28 30  --legend S59-2 S43-2 -o tb201902_RMSfc_comparison.png --rebin 8	
	#python python/plot_producer.py $i$calib_extension -n RMS_fc_0_127_k -b 0  -d 'hist e same' --refuse 28 30  --legend S59-2 S43-2 -o tb201902_RMSfc_0_127_comparison.png --rebin 8
	#python python/plot_producer.py $i$calib_extension -n slope_vs_channel_k  -b 0 -d 'hist e same' --refuse 28 30  --legend S59-2 S43-2 -o tb201902_slope_vs_channel_comparison.png	
	#python python/plot_producer.py $i$calib_extension -n slope_vs_right_strip_k  -b 0 -d 'hist e same' --refuse 28 30  --legend S59-2 S43-2 -o tb201902_slope_vs_right_strip_comparison.png	
	
done

##echo $same
##echo $dates




#for k in ${kpix[*]}
#do
	#kpixname="kpix"
	#kpixname+="$k" 
	##if [$k "=" 28 ] || [$k "=" 30]; then
		##echo "teeeest"
	##fi
	#outputname="$kpixname""_tb201902_pedestal_time_evolv.png"
	#python python/plot_producer.py $same -n pedestals_fc_k$k -b 0 -d 'hist e same' --legend $dates -o $pedestalname --rebin 8
	#outputname="$kpixname""_tb201902_slope_time_evolv.png"
	#python python/plot_producer.py $same -n slope_k$k -b 0  -d 'hist e same'  --legend $dates -o $outputname
	#outputname="$kpixname""_tb201902_RMSfc_time_evolv.png"
	#python python/plot_producer.py  $same -n RMS_fc_k -k $k -b 0  -d 'hist e same'  --legend $dates -o $outputname --rebin 8
#done



outputname="$kpixname""_tb201902_pedestal_time_evolv_new.png"
python python/plot_producer.py /scratch/data/testbeam201902/2019_02_05_16_57_09.bin.newCalib.root $same -n pedestals_fc_k -k 17 22 -b 0 -d 'hist e same' --legend 2019_02_05_16_57_09 $dates -o $pedestalname --rebin 8
outputname="$kpixname""_tb201902_slope_time_evolv_new.png"
python python/plot_producer.py /scratch/data/testbeam201902/2019_02_05_16_57_09.bin.newCalib.root $same -n slope_k -k 17 22 -b 0  -d 'hist e same'  --legend 2019_02_05_16_57_09 $dates -o $outputname
outputname="$kpixname""_tb201902_RMSfc_time_evolv_new.png"
python python/plot_producer.py /scratch/data/testbeam201902/2019_02_05_16_57_09.bin.newCalib.root $same -n RMS_fc_k -k 17 22 -b 0  -d 'hist e same'  --legend 2019_02_05_16_57_09 $dates -o $outputname --rebin 8


