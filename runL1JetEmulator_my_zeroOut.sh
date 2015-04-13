#!/bin/bash
############################################################
## print the command
echo "command called =" $0 $1 $2
## start time of the script
date1=$(date +"%s")         # %s     seconds since 1970-01-01 00:00:00 UTC
#date1Human=$(date -d @$date1)
date1Human=$(date -d @$date1 '+%d/%m/%Y %H:%M:%S')
echo "started on     =" $date1Human
echo "############################################################"
############################################################

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

############################################################
echo "############################################################"
echo "##### log of the SCRIPT =" $0

## finish time and run time of the script
date2=$(date +"%s")         # %s     seconds since 1970-01-01 00:00:00 UTC
#date2Human=$(date -d @$date2)
date2Human=$(date -d @$date2 '+%d/%m/%Y %H:%M:%S')
diff=$(($date2-$date1))
diffHuman=$(date -d @$diff)
#echo $date1
echo "started  on    =" $date1Human
echo "finished on    =" $date2Human
echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo $diffHuman

echo "##### END #####"
############################################################
