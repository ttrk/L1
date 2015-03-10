#!/bin/sh
# compile the macro with g++
g++ L1JetEmulator.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o L1JetEmulator.exe || exit 1

./L1JetEmulator.exe "/export/d00/scratch/luck/Hydjet1p8_2760GeV_L1UpgradeAnalyzer_GT_run1_mc_HIon_L1UpgradeAnalyzer.root" "Hydjet276_JetResults.root" || exit 1
./L1JetEmulator.exe "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_L1UpgradeAnalyzer_GT_MCHI2_74_V3.root" "Hydjet502_JetResults.root" || exit 1

g++ makeRateCurve.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeRateCurve.exe || exit 1

./makeRateCurve.exe "Hydjet276_JetResults.root" || exit 1
./makeRateCurve.exe "Hydjet502_JetResults.root" || exit 1
