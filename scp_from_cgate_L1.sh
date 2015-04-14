#!/bin/sh

## NOT WORKING at the moment
## WARNING : cannot copy directly from scratch to local machine
#http://stackoverflow.com/questions/3598664/creating-a-shell-script-to-run-java-program?rq=1
remoteMachine="tatar@cgate.mit.edu";
sourceFolder="out_L1EmulatorMacros_v3";
sourceDir="output/"$sourceFolder;
sourceFiles="
${sourceDir}/*.pdf
"; # do not copy all sort of files

#source=$remoteMachine":"$sourceFiles;
sourceDir2=$remoteMachine":"$sourceDir;

destination="/home/kaya/Documents/cgate/output/"$sourceFolder;

mkdir -p $destination
# scp "${source}" $destination           # NOT working
scp $sourceDir2"/*.pdf" $destination
