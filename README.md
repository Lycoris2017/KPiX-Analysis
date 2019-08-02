# KPiX-Analysis

Desciption: analysis code for R&D use on the LYCORIS telescope project.

### Notice:
--------------
* in .gitignore: in order to ignore the symlink files, please do
```
find . -type l | sed -e s'/^\.\///g' >> .gitignore
```

### Prerequisite:
--------------

In order to run this package, you will need:
* libkpix.so produced from the kpix DAQ: https://github.com/Lycoris2017/KPiX-Lycoris.git
* you need to either configure the LD_LIBRARY_PATH to have the 'make install'ed kpix package to be used for this package; or you can
```
mkdir lib
cp $KPIX-DAQ/libkpix.so lib/
```

### How-to Install?
--------------

#### Add the central KPiX-Daq repository to get the master branch:

```
git clone https://github.com/Lycoris2017/KPiX-Analysis.git kpix_ana
cd kpix_ana
```

#### (Optional) As developper:

If you are a developper of this software, please click the 'Fork' on top right of this page. Then add your mirror and push your 'dev' branch to it:

```
git remote add myana git@github.com:$YOUR_GITHUB_REPOSITORY/KPiX-Analysis.git
git checkout -b local.dev
git commit -m "init commit, copy from central master branch"
git push -u myana local.dev
```

__VERY IMPORTANT!!!__

When you are confident about your commit, AND you want to MERGE your branch to the central branch: do the following:

1. if you want to avoid spaghetti-history you can __rebase__:
```
git checkout YOUR_COOL_DEV_BRANCH
git rebase master
```
2. then you switch to master, and merge with your development:
```
git checkout master
git merge YOUR_COOL_DEV_BRANCH
```


Collaboration ATTENTION! Please always push to your own branch, then to ask for a pull-request on the central git repository!

#### Compile with CMake

```
mkdir build/; cd build/
cmake ..
make install
```

### How-to Use?
--------------

#### Usage:

If you run an analysis on Acquire binary:
```
bin/analysis acquire.bin
```

If you run ananalysis on Calibration binary:
```
bin/new_CalibrationFitter calib.xml calib.bin
```

If you want the fancy ecal map:
```
cd python
python ecal_mapping.py my_ecal_acquire.bin.root -k 30
```

If you want the coloful strips :
```
cd python
python strip_mapping.py my_strip.bin.root
```

If you want to produce the so-called DA map (the 1024 channels in a 32*(4*8) matrix to activate/deactivate their internal trigger), based on the slope value, pedestal RMS and even hot channels during acquire noise run:
```
python python/channel_disabler.py --calib calib.bin.root --acquire acquire.bin.root -k 26 28
```

or for ecal:
```
python python/channel_disabler.py --calib calib.bin.root --acquire acquire.bin.root -k 30 --ecal True
```

### Repo management
--------------

#### Branch descriptions:

* master: up-to-date latest development version.
* experimental: collaborative branch for development.

#### Subdir structure:

* CMakeLists.txt: cmake config file.
* init.sh: environmental var definition.
* src/ :
```
src/
├── analysis.cxx 
├── analysisEcal.cxx
├── analysisLite.cxx
├── analysis_test.cxx
├── calib_compare.cxx
├── calibrationFitter.cxx
├── count.cxx
├── crossCalib.cpp
├── DuoKpix_CalibFitter.cxx
├── new_CalibrationFitter.cxx
├── ntupleTest.cxx -> ../test/ntupleTest.cxx
├── printDat.cxx
└── subtr_pedestal.cpp
```

* include/ :
```
include/
├── analysis.h
├── ECAL_new1_ng.h
├── kpix_left_and_right.h
├── kpix_left.h
├── kpixmap.h
├── kpix_right.h
├── LinkDef.h
├── ntupleMaker.h
├── RootCompare.h
├── test.h
└── tracker_59_calibration.h
```
* extern/ : external libraries comes here (basically kpix bin decode binaries)
```
extern/
├── include
│   ├── Data.h
│   ├── DataRead.h
│   ├── KpixEvent.h
│   ├── KpixSample.h
│   └── XmlVariables.h
└── src
    ├── Data.cpp
    ├── DataRead.cpp
    ├── KpixEvent.cpp
    ├── KpixSample.cpp
    └── XmlVariables.cpp
```
* scripts/ :
```
scripts/
├── debugger.C
├── mapRun.sh
├── obsolete
│   ├── Makefile.20180830
│   └── test.sh
├── plot_script.sh
├── run.sh
├── tb_analyse_script.sh
└── testbeam_time_diff_script.sh
```

* data/: data needs for `scripts/` or `include/` or `src/`.
```
data/
├── compare.conf
├── kpix_to_tracker.dat
├── tracker_to_kpix_left.dat
└── tracker_to_kpix_right.dat
```

* bin : executables to run analyses or build up tree/particle evnts.
```
bin/
├── analysis
├── analysisLite
├── count
├── DuoKpix_CalibFitter
├── new_CalibrationFitter
├── ntupleTest
└── printDat
```

* test/ : test dir, you do not care.
```
test/
├── ntupleTest.cxx
└── run.sh
```

