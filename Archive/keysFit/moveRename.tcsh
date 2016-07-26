#!/usr/local/bin/tcsh

if ( $# < 2 ) then
    echo The sintaxis is: moveRename.tcsh base new 
    exit -1
endif

set base = $1
set new = $2
set start = `ll ${new}* | wc -l`

foreach file (${base}*.root )
    set bname = `basename $file`
    set newfile = ${new}${start}.root
    mv -v $file $newfile
    @ start++
end


echo Moved
