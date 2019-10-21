# Simple Instructions

### PrintSlopeData.py
--------------

* Description: Print out slopes data from the calibration root file, by searching for histograms called:
`slope_vs_channel_k#_b0`, where # stands for which kpix.


### channel_disabler.py
--------------

* Description: python code to give you a list of 'enabled' or 'disabled' channel list for a kpix
* Usage:
```
python channel_disabler.py [-h] file_in [file_in ...]
```
* Example
```
python channel_disabler.py /dir/2018_06_08_14_31_22.bin.newCalib.root /dir/2018_07_12_19_16_50.bin.root 
```
where the yyyy.root is a pedestal file
* Notice:
  * The __pedestal noisy channels__ is selected based on the normalized events recorded per cycle per channel __< 0.02__.
 

### strip_mapping.py
--------------

* Description: python code to visualize the strip pedestal level or signal response.
* Usage:
```
python strip_mapping.py /dir/2018_07_05_20_02_46.bin.root
```
