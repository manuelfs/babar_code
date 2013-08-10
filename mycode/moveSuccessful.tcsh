#!/bin/tcsh

if ( $# < 2 ) then
    echo The sintaxis is: moveSuccessful.tcsh logdir tcldir
    exit -1
endif

set logdir = $1
set tcldir = $2
set bname = `basename $tcldir`
set outdir = tcls/Done/$bname
mkdir -p $outdir

grep "Successfully" $logdir/* | sed s/"\log:Successfully completed."/tcl/ | awk '{print "basename "$1}' | csh | awk '{print "mv ""'"$tcldir"'"$1" ""'"$outdir"'"}' | csh
#grep "exit code 64" $logdir/* | sed s/"\log:Exited with exit code 64."/tcl/ | awk '{print "basename "$1}' | csh | awk '{print "mv ""'"$tcldir"'"$1" ""'"$outdir"'"}' | csh
#grep "exit code 127" $logdir/* | sed s/"\log:Exited with exit code 127."/tcl/ | awk '{print "basename "$1}' | csh | awk '{print "mv ""'"$tcldir"'"$1" ""'"$outdir"'"}' | csh
#grep "Framework is exiting" $logdir/* | sed s/"\log:Framework is exiting now."/tcl/ | awk '{print "basename "$1}' | csh | awk '{print "mv ""'"$tcldir"'"$1" ""'"$outdir"'"}' | csh

echo Successful files moved to $outdir
