#!/bin/bash

set -x

  # enum seedObject {
  #   emcands,
  #   regions,
  #   subRegions,
  #   twoByTwoJets,
  #   threeByThreeJets
  # };


InputHiForest=("/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_rctconfigNoCuts_HiForestAndEmulatorAndHLT_v7.root" "/mnt/hadoop/cms/store/user/luck/L1Emulator/AllQCDPhoton30_PhotonFilter20GeV_eta3_TuneZ2_PbPb_5020GeV_embedded_rctconfigNoCuts_HiForest.root")

SampleType=(Hydjet502 Photon502)

SeedTypes=(emcands regions subregions 2x2jets 3x3jets)

EtaCuts=(1.44 2.0)

OfflineIsolation=(0 1)

# compile the macros with g++
g++ makeTurnOn_fromSameFile_photons.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn_fromSameFile_photons.exe || exit 1

for SAMPLE in 0 1
do
    for ALGO in 0 1 2 3 4
    do
	for ETA in 0 #1
	do
	    for ISO in 0 1
	    do
		isolabel=""
		if [[ $ISO == 1 ]]
		then
		    isolabel="iso"
		fi
		output="hist_${SampleType[SAMPLE]}_${isolabel}photons_eta${EtaCuts[ETA]}_${SeedTypes[ALGO]}.root"
		./makeTurnOn_fromSameFile_photons.exe "${InputHiForest[SAMPLE]}" "$output" "$ALGO" "${EtaCuts[ETA]}" "$ISO"
	    done
	done
    done
done
