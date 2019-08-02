#------------------
#- Author: Mengqing Wu <mengqing.wu@desy.de>
#- October 1st, 2018
#- for LYCORIS project at DESY.
#------------------

#!bin/sh

# to be fixed: need to manually add the Linkdef files for RPATH (typedef or classes for TTree)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib

# want to access executable without type its dir?
export PATH=$PATH:$PWD/bin
