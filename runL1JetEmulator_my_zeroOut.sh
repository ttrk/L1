#!/bin/bash

# The input files are listed here : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TValidationSamples

# # Alex's email that led to this bash script
# # Kaya, Bi Ran,

# # I've uploaded a few macros that you can use to draw turn-on curves based on different L1 seeds for offline photons:

# # For use with most datasets:
# # https://github.com/CmsHI/L1EmulatorMacros/blob/master/makeTurnOn_photons.C
# # You should find the L1Tuples/HiForests used as input here:
# # https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TValidationSamples

# # For use only with the 5.02 TeV Hydjet sample:
# # https://github.com/CmsHI/L1EmulatorMacros/blob/master/makeTurnOn_fromSameFile_photons.C
# # ...


# "/mnt/hadoop/cms/store/user/luck/L1Emulator/HiForest_PbPb_photon2030.root"
# DATA   = Skims
# sample = Photon20/30

# "/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_HiForestAndEmulator_v3.root"
# DATA   = MC
# sample = 5.02TeV MB Hydjet
InputHiForest=(
"/mnt/hadoop/cms/store/user/luck/L1Emulator/HiForest_PbPb_photon2030.root"
"/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_HiForestAndEmulator_v3.root"
);
InputL1=(
"/export/d00/scratch/ginnocen/PhotonSamples/HIHighPt-HIRun2011-RAW-photon20and30_v2_GR_P_V27A_L1UpgradeAnalyzer.root"
);

outDirectory="/export/d00/scratch/tatar/output/out_L1EmulatorMacros_v3";

outFileNameHist=(
"hist_Photon2030_zeroOut_NxN.root"
"hist_HydjetMB_502TeV_zeroOut_NxN.root"
);
outFilePathHist=(
$outDirectory"/"${outFileNameHist[0]}
$outDirectory"/"${outFileNameHist[1]}
);

outFileNameHist2x2=(
"hist_Photon2030_zeroOut_2x2.root"
"hist_HydjetMB_502TeV_zeroOut_2x2.root"
);
outFilePathHist2x2=(
$outDirectory"/"${outFileNameHist2x2[0]}
$outDirectory"/"${outFileNameHist2x2[1]}
);

outFileNameHist3x3=(
"hist_Photon2030_zeroOut_3x3.root"
"hist_HydjetMB_502TeV_zeroOut_3x3.root"
);
outFilePathHist3x3=(
$outDirectory"/"${outFileNameHist3x3[0]}
$outDirectory"/"${outFileNameHist3x3[1]}
);

# compile the macros with g++
g++ makeTurnOn_photons.C              $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn_photons.exe              || exit 1
g++ makeTurnOn_fromSameFile_photons.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn_fromSameFile_photons.exe || exit 1
##g++ plotTurnOn.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o plotTurnOn.exe || exit 1

for sampleNum in 0 1
do
	if [ $sampleNum -eq 0 ]
	then
	    ./makeTurnOn_photons.exe "${InputL1[sampleNum]}" "${InputHiForest[sampleNum]}" "${outFilePathHist2x2[sampleNum]}" || exit 1
	elif [ $sampleNum -eq 1 ]
	then
	    ./makeTurnOn_fromSameFile_photons.exe            "${InputHiForest[sampleNum]}" "${outFilePathHist2x2[sampleNum]}" || exit 1
        fi
done

for sampleNum in 0 1
do
	if [ $sampleNum -eq 0 ]
	then
	    ./makeTurnOn_photons.exe "${InputL1[sampleNum]}" "${InputHiForest[sampleNum]}" "${outFilePathHist3x3[sampleNum]}" || exit 1
	elif [ $sampleNum -eq 1 ]
	then
	    ./makeTurnOn_fromSameFile_photons.exe            "${InputHiForest[sampleNum]}" "${outFilePathHist3x3[sampleNum]}" || exit 1
        fi
done














# # set -x

# # This is a master file to run the emulator macro and then print out rates on multiple datasets at once
# # You can find the up-to-date list of samples here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TValidationSamples


# for sampleNum in 0 1 2
# do
#     for algo in 7 8 9 10 #0 1 2 3 4 5 6 7 8 9 10 # only run new combinations now
#     do
# 	L1Output="~/scratch/EmulatorResults/${InputType[sampleNum]}_JetResults_${AlgoVariations[algo]}.root"
# 	HistOutput="hist_${InputType[sampleNum]}_${AlgoVariations[algo]}.root"
# 	PlotOutputTag="${InputType[sampleNum]}_${AlgoVariations[algo]}"
# 	./L1JetEmulator.exe "${InputL1[sampleNum]}" "$L1Output" $algo || exit 1
# 	./makeRateCurve.exe "$L1Output" 1 || exit 1
# 	if [[ $sampleNum -eq 0 ]] || [[ $sampleNum -eq 1 ]]
# 	then
# 	   ./makeTurnOn.exe "$L1Output" "${InputHiForest[sampleNum]}" "$HistOutput" || exit 1
# 	elif [[ $sampleNum -eq 2 ]]
# 	then
# 	   ./makeTurnOn_fromSameFile.exe "$L1Output" "${InputHiForest[sampleNum]}" "$HistOutput" || exit 1
# 	fi
# 	./plotTurnOn.exe "$HistOutput" "$PlotOutputTag" || exit 1
#     done
# done

