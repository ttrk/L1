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

outDirectory="/export/d00/scratch/tatar/output/out_L1EmulatorMacros_v3";

outFileNameHist2x2=(
"hist_Photon2030_zeroOut_2x2"
"hist_HydjetMB_502TeV_zeroOut_2x2"
);
outFilePathHist2x2=(
$outDirectory"/"${outFileNameHist2x2[0]}".root"
$outDirectory"/"${outFileNameHist2x2[1]}".root"
);
outFilePathTag2x2=(
$outDirectory"/"${outFileNameHist2x2[0]}
$outDirectory"/"${outFileNameHist2x2[1]}
);

outFileNameHist3x3=(
"hist_Photon2030_zeroOut_3x3"
"hist_HydjetMB_502TeV_zeroOut_3x3"
);
outFilePathHist3x3=(
$outDirectory"/"${outFileNameHist3x3[0]}".root"
$outDirectory"/"${outFileNameHist3x3[1]}".root"
);
outFilePathTag3x3=(
$outDirectory"/"${outFileNameHist3x3[0]}
$outDirectory"/"${outFileNameHist3x3[1]}
);

# compile the macros with g++
g++ plotTurnOn_photons.C   $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o plotTurnOn_photons.exe     || exit 1
##g++ plotTurnOn.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o plotTurnOn.exe || exit 1

for sampleNum in 0 1
do
    ./plotTurnOn_photons.exe "${outFilePathHist2x2[sampleNum]}" "${outFilePathTag2x2[sampleNum]}"  || exit 1
done

for sampleNum in 0 1
do
    ./plotTurnOn_photons.exe "${outFilePathHist3x3[sampleNum]}" "${outFilePathTag3x3[sampleNum]}" || exit 1
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

