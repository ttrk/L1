#!/bin/sh

#http://stackoverflow.com/questions/3598664/creating-a-shell-script-to-run-java-program?rq=1

remoteMachine="tatar@cgate.mit.edu";
sourceFolder="out_L1EmulatorMacros";
sourceDir="/export/d00/scratch/tatar/output/"$sourceFolder;
sourceFiles="
${sourceDir}/*.gif
";         # do not copy all sort of files
source=$remoteMachine":"$sourceFiles;

destination="~/Documents/cgate/output/"$sourceFolder;

mkdir $destination
scp $sourceFiles $destination
