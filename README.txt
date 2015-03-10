Set of macros to run L1 Jet algorithm on top of previous L1 results.

Instructions:
The macro L1JetEmulator.C contains a copy of the L1 jet emulator. It expects as input a file made by the L1UpgradeAnalyzer code, but only uses the region information. To test a change to the L1 jet algorithm, modify this macro and use runL1JetEmulator.sh to run it over both the Hydjet 2.76 and 5.02 samples.

Once the output is created, you can use the makeRateCurve.C macro to draw the L1 Jet rate as a function of physical pT in GeV.

The macro drawTurnOnCurves.C will also draw a set of Turn-on curves and requires event matching between a HiForest file and the output of the L1 Jet emulator. You can change the L1 thresholds inside the macro.
