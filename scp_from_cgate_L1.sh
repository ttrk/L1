#!/bin/sh

#http://stackoverflow.com/questions/3598664/creating-a-shell-script-to-run-java-program?rq=1

remoteMachine="tatar@cgate.mit.edu";
sourceDir="~/output/out_L1EmulatorMacros";
sourceFiles="
${sourceDir}/*.C
${sourceDir}/*.h
${sourceDir}/*.sh
${sourceDir}/*.txt
";         # do not copy all sort of files
source=$remoteMachine":"$sourceFiles;

destination="~/Documents/cgate/output/";

scp $sourceFiles $destination
