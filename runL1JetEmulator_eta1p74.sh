#!/bin/sh

# This is a master file to run the emulator macro and then print out rates on multiple datasets at once
# You can find the up-to-date list of samples here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TValidationSamples

# Use this to make sure that you do not mix up different algorithm tweaks.
#TAGNAME="" #for nominal running
#TAGNAME="_zeroWalls" # for zeroing out walls
#TAGNAME="_doubleIteration"
#TAGNAME="_1x1Jets"
#TAGNAME="_truncatedMean"
#TAGNAME="_etaReflection"
TAGNAME="_eta1p74"

# compile the macro with g++
g++ L1JetEmulator${TAGNAME}.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o L1JetEmulator${TAGNAME}.exe || exit 1

./L1JetEmulator${TAGNAME}.exe "/export/d00/scratch/luck/Hydjet1p8_2760GeV_L1UpgradeAnalyzer_GT_run1_mc_HIon_L1UpgradeAnalyzer.root" "~/scratch/Hydjet276_JetResults${TAGNAME}.root" || exit 1
./L1JetEmulator${TAGNAME}.exe "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_L1UpgradeAnalyzer_GT_MCHI2_74_V3.root" "~/scratch/Hydjet502_JetResults${TAGNAME}.root" || exit 1
./L1JetEmulator${TAGNAME}.exe "/mnt/hadoop/cms/store/user/men1/L1Data/HIL1DPG/MinBias/HIMinBiasUPC_Skim_HLT_HIMinBiasHfOrBSC_v2_CaloRegionEta516_CMSSW740pre7/L1NTupleMBHI.root" "~/scratch/MBData_JetResults${TAGNAME}.root" || exit 1

g++ makeRateCurve.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeRateCurve.exe || exit 1

./makeRateCurve.exe "~/scratch/Hydjet276_JetResults${TAGNAME}.root" || exit 1
./makeRateCurve.exe "~/scratch/Hydjet502_JetResults${TAGNAME}.root" || exit 1
./makeRateCurve.exe "~/scratch/MBData_JetResults${TAGNAME}.root" || exit 1

g++ makeTurnOn${TAGNAME}.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn${TAGNAME}.exe || exit 1

./makeTurnOn${TAGNAME}.exe "~/scratch/Hydjet276_JetResults${TAGNAME}.root" "/mnt/hadoop/cms/store/user/ginnocen/Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV/HiMinBias_Forest_26June2014/d9ab4aca1923b3220eacf8ee0d550950/*.root" "hist_hydjet276${TAGNAME}.root" 1 || exit 1

./makeTurnOn${TAGNAME}.exe "~/scratch/Hydjet502_JetResults${TAGNAME}.root" "/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_740pre8_MCHI2_74_V3_53XBS_HiForest.root" "hist_hydjet502${TAGNAME}.root" 1 || exit 1

./makeTurnOn${TAGNAME}.exe "~/scratch/MBData_JetResults${TAGNAME}.root" "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged_v2/HiForest_PbPb_Data_minbias_fromSkim_v3.root" "hist_MBData${TAGNAME}.root" 0 || exit 1


g++ plotTurnOn.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o plotTurnOn.exe || exit 1

./plotTurnOn.exe "hist_hydjet276${TAGNAME}.root" "hydjet276${TAGNAME}" || exit 1
./plotTurnOn.exe "hist_hydjet502${TAGNAME}.root" "hydjet502${TAGNAME}" || exit 1
./plotTurnOn.exe "hist_MBData${TAGNAME}.root" "MBData${TAGNAME}" || exit 1
