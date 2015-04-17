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

# http://stackoverflow.com/questions/8880603/loop-through-array-of-strings-in-bash-script
## declare an array variable for the types of algorithms
declare -a algoTypes=("NOzeroOut_BkgSubtract" "zeroOut_2x2" "zeroOut_3x3");

fileNameHist=(
"hist_Photon2030"
"hist_HydjetMB_502TeV"
);

# compile the macros with g++
g++ plotTurnOn_photons.C   $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o plotTurnOn_photons.exe     || exit 1

for sampleNum in 0 1
do
    for i in "${algoTypes[@]}"
    do
	filePathHist=$outDirectory"/"${fileNameHist[sampleNum]}"_"${i}".root"
	outFilePathTag=$outDirectory"/"${fileNameHist[sampleNum]}"_"${i}
	./plotTurnOn_photons.exe "${filePathHist}" "${outFilePathTag}"  || exit 1
    done
done

