#!/bin/sh

# This is a master file to run the emulator macro and then print out rates on multiple datasets at once
# You can find the up-to-date list of samples here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TValidationSamples

# compile the macro with g++
g++ L1JetEmulator.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o L1JetEmulator.exe || exit 1

### Kayas modification
outDirectory="/net/hisrv0001/home/tatar/output/out_L1EmulatorMacros";
outFileName1="Hydjet276_JetResults.root";
outFileName2="Hydjet502_JetResults.root";
outFileName1Hist="hist_hydjet276.root";
### Kayas modification

./L1JetEmulator.exe "/export/d00/scratch/luck/Hydjet1p8_2760GeV_L1UpgradeAnalyzer_GT_run1_mc_HIon_L1UpgradeAnalyzer.root" $outDirectory"/"$outFileName1 || exit 1
./L1JetEmulator.exe "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_L1UpgradeAnalyzer_GT_MCHI2_74_V3.root" $outDirectory"/"$outFileName2 || exit 1

g++ makeRateCurve.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeRateCurve.exe || exit 1

./makeRateCurve.exe $outDirectory"/"$outFileName1 || exit 1
./makeRateCurve.exe $outDirectory"/"$outFileName2 || exit 1

g++ makeTurnOn.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn.exe || exit 1

./makeTurnOn.exe $outDirectory"/"$outFileName1 "/mnt/hadoop/cms/store/user/ginnocen/Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV/HiMinBias_Forest_26June2014/d9ab4aca1923b3220eacf8ee0d550950/*.root" $outDirectory"/"$outFileName1Hist  1 || exit 1
