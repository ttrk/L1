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
"/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_HiForestAndEmulator_v5.root"
"/mnt/hadoop/cms/store/user/luck/L1Emulator/AllQCDPhoton30_PhotonFilter20GeV_eta3_TuneZ2_PbPb_5020GeV_actuallyEmbedded_HiForest.root"
"/mnt/hadoop/cms/store/user/luck/L1Emulator/PyquenUnquenched_DiJet_pt30_PbPb_5020GeV_actuallyEmbedded_HiForest.root"
);

outDirectory="/export/d00/scratch/tatar/output/out_L1EmulatorMacros_v3";

declare -a algoTypes=("NOzeroOut_BkgSubtract" "zeroOut_2x2" "zeroOut_3x3");
#declare -a algoTypes=("zeroOut_1x1");

outFileNamePrefix=(
"hist_HydjetMB_502TeV_v5"
"hist_AllQCDPhoton30_PhotonFilter20GeV"
"hist_PyquenUnquenched_DiJet_pt30_PbPb_5020GeV"
);

# compile the macros with g++
g++ makeTurnOn_fromSameFile_photons.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn_fromSameFile_photons.exe || exit 1

for sampleNum in 0 1 2
do
    for i in "${algoTypes[@]}"
    do
	outFileNameHist=${outFileNamePrefix[sampleNum]}"_"${i}".root"
	outFilePathHist=$outDirectory"/"$outFileNameHist
	./makeTurnOn_fromSameFile_photons.exe            "${InputHiForest[sampleNum]}" "${outFilePathHist}" || exit 1
    done
done
