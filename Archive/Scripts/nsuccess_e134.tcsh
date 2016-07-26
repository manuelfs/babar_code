#!/bin/tcsh

    if ( $# < 1 ) then
	echo I need the folder
	exit -1
    endif

set folder = $1
set StemFolder = `basename $folder`
set d = `date | awk '{print $6 $2 $3}'`
printf "\n"
foreach dir ( $folder/* )
    set bname = `basename $dir`
    set nFiles = `ls tcls/${StemFolder}/${bname}/ | wc -l`
    set e64 = `find ${dir}/ -name "*.log" | xargs grep "exit code 64" | wc -l`
    set e134 = `find ${dir}/ -name "*.log" | xargs grep "exit code 134" | wc -l`
    set n = `find ${dir}/ -name "*.log" | xargs grep "Successfully completed" | wc -l`
    set Exiting = `find ${dir}/ -name "*.log" | xargs grep "Framework is" | wc -l`
    set f = `find ${dir}/ -name "*.log" | xargs grep "Exited" | wc -l`
    set f = `expr $f / 2`
    echo ${bname} - Total ${nFiles}/${Exiting}: $n successes, ${f}/${e134} exits/e134 - Exit 64: $e64
end
printf "\n"
