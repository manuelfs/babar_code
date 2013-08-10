#!/bin/tcsh

if ( $# < 1 ) then
    echo The sintaxis is: moveBack.tcsh AWG82/logfiles/ntuples/R24/SP-1...
    exit -1
endif


set folder = $1
set bname = `basename $folder`
set dir = `dirname $folder`
set StemFolder = `basename $dir`
set dir = `dirname $dir`
set ntuples = AWG82/ntuples/$bname

mv $dir/$bname/* $folder/
rm -r $dir/$bname
mv $ntuples/* AWG82/ntuples/$StemFolder/$bname/
rm -r $ntuples
mv tcls/$StemFolder/$bname/* tcls/Done/$bname/
rm -r tcls/$StemFolder/$bname/ 
mv tcls/Done/$bname tcls/$StemFolder/









