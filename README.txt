Set of macros to run L1 Jet algorithm on top of previous L1 results.

Instructions:
The macro L1JetEmulator.C contains a copy of the L1 jet emulator. It expects as input a file made by the L1UpgradeAnalyzer code, but only uses the region information. To test a change to the L1 jet algorithm, modify this macro and use runL1JetEmulator.sh to run the algorithm, make rate curves, and make turn-on curves.

makeRateCurve.C macro will draw the L1 Jet rate as a function of l1 physical pT in GeV. It can be run interactively in a root session to draw the rate curve, or in batch mode which will dump the rates to text.

The macro makeTurnOn.C will match events from the algorithm output and a HiForest and then make a root file with turn on curves in it. You can use the macro plotTurnOn.C interactively in a root session to then make pretty plots out of the results.
