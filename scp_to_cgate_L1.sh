#!/bin/sh

#http://stackoverflow.com/questions/3598664/creating-a-shell-script-to-run-java-program?rq=1

#sourceDir="L1EmulatorMacros"
# sourceFiles="
# *.C
# *.h
# *.sh
# *.txt
# ";         # do not copy all sort of files
sourceFiles="makeTurnOn_vIso_Xjg.C";
remoteMachine="tatar@cgate.mit.edu";
destinationInRemote="~/code";
destination=$remoteMachine":"$destinationInRemote;

scp $sourceFiles $destination
