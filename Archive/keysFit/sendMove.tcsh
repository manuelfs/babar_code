#!/usr/local/bin/tcsh

set i = 0
while ( $i <= 19 )
    foreach smoo (70 100 130)
	set rot =  AWG82/results/keys/Archive/root/new/hKeysRot${i}_${smoo}_
	set new =  AWG82/results/keys/Archive/root/new/hKeysNoR${i}_${smoo}_
	# babar_code/keysFit/moveRename.tcsh $base $new
	set start = `ll ${rot}* | wc -l`
	echo Rot ${i}_${smoo} has $start
	set start = `ll ${new}* | wc -l`
	echo noRot ${i}_${smoo} has $start
    end
    @ i = $i + 1
end
