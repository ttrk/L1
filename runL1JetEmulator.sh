#!/bin/sh

# This is a master file to run the emulator macro and then print out rates on multiple datasets at once
# You can find the up-to-date list of samples here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TValidationSamples

# Use this to make sure that you do not mix up different algorithm tweaks.
#TAGNAME=""                  algo=""               #for nominal running
#TAGNAME="_zeroWalls"        algo="zeroWalls"      # for zeroing out walls
#TAGNAME="_doubleIteration"  algo="doubleIteration";
#TAGNAME="_smallJets"        algo="smallJets";
#TAGNAME="_truncatedMean"    algo="truncatedMean";
TAGNAME="_etaReflection";    algo="etaReflection";

### Kayas modification
outDirectory="/export/d00/scratch/tatar/output/out_L1EmulatorMacros_v2"; 	# do not end this with "/"

outFileName1="Hydjet276_JetResults${TAGNAME}.root";
outFileName2="Hydjet502_JetResults${TAGNAME}.root";
outFileName3="MBData_JetResults${TAGNAME}.root";

prefix1="hydjet276";
prefix2="hydjet502";
prefix3="MBData";

outFileName1Hist="hist_${prefix1}${TAGNAME}.root";
outFileName2Hist="hist_${prefix2}${TAGNAME}.root";
outFileName3Hist="hist_${prefix3}${TAGNAME}.root";

outFilePath1=$outDirectory"/"$outFileName1;
outFilePath2=$outDirectory"/"$outFileName2;
outFilePath3=$outDirectory"/"$outFileName3;

outFilePath1Hist=$outDirectory"/"$outFileName1Hist;
outFilePath2Hist=$outDirectory"/"$outFileName2Hist;
outFilePath3Hist=$outDirectory"/"$outFileName3Hist;
### Kayas modification

# compile the macro with g++
g++ L1JetEmulator.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o L1JetEmulator.exe || exit 1

./L1JetEmulator.exe "/export/d00/scratch/luck/Hydjet1p8_2760GeV_L1UpgradeAnalyzer_GT_run1_mc_HIon_L1UpgradeAnalyzer.root" $outFilePath1 $algo || exit 1
./L1JetEmulator.exe "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_L1UpgradeAnalyzer_GT_MCHI2_74_V3.root" $outFilePath2 $algo || exit 1
./L1JetEmulator.exe "/mnt/hadoop/cms/store/user/men1/L1Data/HIL1DPG/MinBias/HIMinBiasUPC_Skim_HLT_HIMinBiasHfOrBSC_v2_CaloRegionEta516_CMSSW740pre7/L1NTupleMBHI.root" $outFilePath3 $algo || exit 1

g++ makeRateCurve.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeRateCurve.exe || exit 1

./makeRateCurve.exe $outFilePath1 > "${outDirectory}/rates30kHz_${prefix1}${TAGNAME}.txt"   || exit 1
./makeRateCurve.exe $outFilePath2 > "${outDirectory}/rates30kHz_${prefix2}${TAGNAME}.txt"  || exit 1
./makeRateCurve.exe $outFilePath3 > "${outDirectory}/rates30kHz_${prefix3}${TAGNAME}.txt"  || exit 1

g++ makeTurnOn.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn.exe || exit 1

./makeTurnOn.exe $outFilePath1 "/mnt/hadoop/cms/store/user/ginnocen/Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV/HiMinBias_Forest_26June2014/d9ab4aca1923b3220eacf8ee0d550950/*.root" $outFilePath1Hist 1 || exit 1

./makeTurnOn.exe $outFilePath2 "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_HiForest_partial.root" $outFilePath2Hist 1 || exit 1

./makeTurnOn.exe $outFilePath3 "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged_v2/HiForest_PbPb_Data_minbias_fromSkim_v3.root" $outFilePath3Hist 0 || exit 1


g++ plotTurnOn.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o plotTurnOn.exe || exit 1

./plotTurnOn.exe $outFilePath1Hist "${outDirectory}/${prefix1}${TAGNAME}" || exit 1
./plotTurnOn.exe $outFilePath2Hist "${outDirectory}/${prefix2}${TAGNAME}" || exit 1
./plotTurnOn.exe $outFilePath3Hist "${outDirectory}/${prefix3}${TAGNAME}" || exit 1

