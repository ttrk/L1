#!/bin/sh

#http://stackoverflow.com/questions/3598664/creating-a-shell-script-to-run-java-program?rq=1

#sourceDir="L1EmulatorMacros"
sourceFiles="
*.C
*.h
*.sh
*.txt
";         # do not copy all sort of files
remoteMachine="tatar@cgate.mit.edu";
destinationInRemote="~/CMSSW_7_4_0_pre8/src/L1EmulatorMacros";
destination=$remoteMachine":"$destinationInRemote;

scp $sourceFiles $destination
