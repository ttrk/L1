#!/bin/sh

#http://stackoverflow.com/questions/3598664/creating-a-shell-script-to-run-java-program?rq=1
remoteMachine="tatar@cgate.mit.edu";
sourceFolder="out_L1EmulatorMacros";
sourceDir="/export/d00/scratch/tatar/output/"$sourceFolder;
sourceFiles="
${sourceDir}/*.pdf
"; # do not copy all sort of files
source=$remoteMachine":"$sourceFiles;

destination="/home/kaya/Documents/cgate/output/"$sourceFolder;

mkdir $destination
scp $source $destination


###!/bin/sh

###http://stackoverflow.com/questions/3598664/creating-a-shell-script-to-run-java-program?rq=1
##
##remoteMachine="tatar@cgate.mit.edu";
##sourceFolder="out_L1EmulatorMacros";
##sourceDir="/export/d00/scratch/tatar/output/"$sourceFolder;
##sourceFiles="
##${sourceDir}/*.pdf
##";         # do not copy all sort of files
###source=$remoteMachine":"$sourceFiles;
##
##destination2="~/output/"$sourceFolder;
##source2=$remoteMachine":"$destination2;
##
##destination="/home/kaya/Documents/cgate/output/"$sourceFolder;
##
###connect to remote
##ssh -24Y tatar@cgate.mit.edu
### on remote machine
##cp    $sourceFiles $destination2
##exit 1
### on local machine
##mkdir $destination
##scp   $source2     $destination

