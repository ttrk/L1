#!/bin/bash

set -x

# This is a master file to run the emulator macro and then print out rates on multiple datasets at once
# You can find the up-to-date list of samples here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TValidationSamples

  # enum algoVariation {
  #0   nominal,
  #1   zeroWalls,
  #2   doubleSubtraction,
  #3   sigmaSubtraction,
  #4   barrelOnly,
  #5   oneByOne,
  #6   twoByTwo,
  #7  oneByOneANDzeroWalls,
  #8  oneByOneANDzeroWallsANDsigmaSubtraction,
  #9  twoByTwoANDzeroWalls,
  #10  twoByTwoANDzeroWallsANDsigmaSubtraction
  # };

InputType=(MBData Hydjet276 Hydjet502)
InputHiForest=("/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged_v2/HiForest_PbPb_Data_minbias_fromSkim_v3.root" "/mnt/hadoop/cms/store/user/ginnocen/Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV/HiMinBias_Forest_26June2014/d9ab4aca1923b3220eacf8ee0d550950/*.root" "/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_HiForestAndEmulator_v3.root")
InputL1=("/mnt/hadoop/cms/store/user/men1/L1Data/HIL1DPG/MinBias/HIMinBiasUPC_Skim_HLT_HIMinBiasHfOrBSC_v2_CaloRegionEta516_CMSSW740pre7/L1NTupleMBHIFS.root" "/export/d00/scratch/luck/Hydjet1p8_2760GeV_L1UpgradeAnalyzer_GT_run1_mc_HIon_L1UpgradeAnalyzer.root" "/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_HiForestAndEmulator_v3.root")

AlgoVariations=(nominal zeroWalls doubleSubtraction sigmaSubtraction barrelOnly oneByOne twoByTwo oneByOneAndzeroWalls oneByOneANDzeroWallsANDsigmaSubtraction twoByTwoANDzeroWalls twoByTwoANDzeroWallsANDsigmaSubtraction)

# compile the macros with g++
g++ L1JetEmulator.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o L1JetEmulator.exe || exit 1
g++ makeRateCurve.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeRateCurve.exe || exit 1
g++ makeTurnOn.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn.exe || exit 1
g++ makeTurnOn_fromSameFile.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn_fromSameFile.exe || exit 1
g++ plotTurnOn.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o plotTurnOn.exe || exit 1


for sampleNum in 0 1 2
do
    for algo in 7 8 9 10 #0 1 2 3 4 5 6 7 8 9 10 # only run new combinations now
    do
	L1Output="~/scratch/EmulatorResults/${InputType[sampleNum]}_JetResults_${AlgoVariations[algo]}.root"
	HistOutput="hist_${InputType[sampleNum]}_${AlgoVariations[algo]}.root"
	PlotOutputTag="${InputType[sampleNum]}_${AlgoVariations[algo]}"
	./L1JetEmulator.exe "${InputL1[sampleNum]}" "$L1Output" $algo || exit 1
	./makeRateCurve.exe "$L1Output" 1 || exit 1
	if [[ $sampleNum -eq 0 ]] || [[ $sampleNum -eq 1 ]]
	then
	   ./makeTurnOn.exe "$L1Output" "${InputHiForest[sampleNum]}" "$HistOutput" || exit 1
	elif [[ $sampleNum -eq 2 ]]
	then
	   ./makeTurnOn_fromSameFile.exe "$L1Output" "${InputHiForest[sampleNum]}" "$HistOutput" || exit 1
	fi
	#./plotTurnOn.exe "$HistOutput" "$PlotOutputTag" || exit 1
    done
done
