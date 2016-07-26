#!/usr/local/bin/tcsh

set wFFBFSig = babar_code/Reweight/wFFBFSig.txt
set wFFBF = babar_code/Reweight/wFFBF.txt
set q = xlong
set x = xxl

set d = `date +%y-%m-%d`
set logdir = AWG82/logfiles/results/keys/single/$d/
mkdir -p $logdir
bsub -q $q -o $logdir/logkeys_0.log FitPdfSampleKEYS   0  85 1 1000 300 1 0.8  1.95  45 $wFFBF         #Minimum at 0.86
bsub -q $q -o $logdir/logkeys_1.log FitPdfSampleKEYS   1 160 1 1000 300 1 -1   1.9   80 $wFFBF         #Minimum at 0.75
bsub -q $q -o $logdir/logkeys_2.log FitPdfSampleKEYS   2  90 1 1000 300 1 1    1.85  60 $wFFBF         #Minimum at 0.86
bsub -q $q -o $logdir/logkeys_3.log FitPdfSampleKEYS   3 105 1 1000 300 1 0.45 1.5   55 $wFFBF         #Minimum at 0.98
bsub -q $q -o $logdir/logkeys_4.log FitPdfSampleKEYS   4 135 1 1000 300 1 0.4  1.5   60 $wFFBF         #Minimum at 0.98
bsub -q $q -o $logdir/logkeys_5.log FitPdfSampleKEYS   5 130 1 1000 300 1 0.4  1.45  60 $wFFBF         #Minimum at 1.33
bsub -q $q -o $logdir/logkeys_6.log FitPdfSampleKEYS   6 180 1 1000 300 1 0.4  1.75  75 $wFFBFSig      #Minimum at 0.98
bsub -q $q -o $logdir/logkeys_7.log FitPdfSampleKEYS   7 180 1 1000 300 1 0.6  1.7   75 $wFFBFSig      #Minimum at 0.86
bsub -q $q -o $logdir/logkeys_8.log FitPdfSampleKEYS   8 190 1 1000 300 1 0.75 1.9  100 $wFFBF         #Minimum at 1.54
bsub -q $q -o $logdir/logkeys_9.log FitPdfSampleKEYS   9 180 1 1000 300 1 0.4  1.45  85 $wFFBF         #Minimum at 1.31
bsub -q $q -o $logdir/logkeys_10.log FitPdfSampleKEYS 10  95 1 1000 300 1 0.75 1.9   55 $wFFBF         #Minimum at 0.86
bsub -q $q -o $logdir/logkeys_11.log FitPdfSampleKEYS 11 180 1 1000 300 1 -1   1.9   75 $wFFBF         #Minimum at 0.75
bsub -q $q -o $logdir/logkeys_12.log FitPdfSampleKEYS 12 100 1 1000 300 1 0.85 1.9   55 $wFFBF         #Minimum at 0.86
bsub -q $q -o $logdir/logkeys_13.log FitPdfSampleKEYS 13 100 1 1000 300 1 0.4  1.6   55 $wFFBF         #Minimum at 0.97
bsub -q $q -o $logdir/logkeys_14.log FitPdfSampleKEYS 14 120 1 1000 300 1 0.4  1.5   60 $wFFBF         #Minimum at 0.86
bsub -q $q -o $logdir/logkeys_15.log FitPdfSampleKEYS 15 120 1 1000 300 1 0.4  1.5   60 $wFFBF         #Minimum at 1.18
bsub -q $q -o $logdir/logkeys_16.log FitPdfSampleKEYS 16 150 1 1000 300 1 0.5  1.85  70 $wFFBFSig      #Minimum at 0.86
bsub -q $q -o $logdir/logkeys_17.log FitPdfSampleKEYS 17 150 1 1000 300 1 0.5  1.85  70 $wFFBFSig      #Minimum at 0.98
bsub -q $q -o $logdir/logkeys_18.log FitPdfSampleKEYS 18 140 1 1000 300 1 -1   -1   100 $wFFBF         #Minimum at 1.06
bsub -q $q -o $logdir/logkeys_19.log FitPdfSampleKEYS 19 170 1 1000 300 1 0.35 1.45  85 $wFFBF         #Minimum at 1.49
bsub -q $q -o $logdir/logkeys_20.log FitPdfSampleKEYS 20 120 1 1000 300 1 0.8  1.8   70 $wFFBFSig      #Minimum at 0.64
bsub -q $q -o $logdir/logkeys_21.log FitPdfSampleKEYS 21 130 1 1000 300 1 0.5  1.9   80 $wFFBFSig      #Minimum at 1.09
bsub -q $q -o $logdir/logkeys_22.log FitPdfSampleKEYS 22 150 1 1000 300 1 0.5  1.8   85 $wFFBFSig      #Minimum at 0.98
bsub -q $q -o $logdir/logkeys_23.log FitPdfSampleKEYS 23 100 1 1000 300 1 0.7  1.8   70 $wFFBFSig      #Minimum at 0.64
bsub -q $q -o $logdir/logkeys_24.log FitPdfSampleKEYS 24 120 1 1000 300 1 0.5  2     80 $wFFBFSig      #Minimum at 0.95
bsub -q $q -o $logdir/logkeys_25.log FitPdfSampleKEYS 25 100 1 1000 300 1 -1   -1   130 $wFFBFSig      #Minimum at 0.95
bsub -q $q -o $logdir/logkeys_26.log FitPdfSampleKEYS 26 110 1 1000 300 1 0.55 1.8   70 $wFFBFSig      #Minimum at 0.53
bsub -q $q -o $logdir/logkeys_27.log FitPdfSampleKEYS 27 120 1 1000 300 1 0.4  1.95  75 $wFFBFSig      #Minimum at 0.73
bsub -q $q -o $logdir/logkeys_28.log FitPdfSampleKEYS 28 140 1 1000 300 1 0.4  1.9   80 $wFFBFSig      #Minimum at 0.98
bsub -q $q -o $logdir/logkeys_29.log FitPdfSampleKEYS 29 130 1 1000 300 1 0.6  1.85  80 $wFFBFSig      #Minimum at 0.53
bsub -q $q -o $logdir/logkeys_30.log FitPdfSampleKEYS 30 140 1 1000 300 1 0.4  1.9   85 $wFFBFSig      #Minimum at 0.53
bsub -q $q -o $logdir/logkeys_31.log FitPdfSampleKEYS 31 150 1 1000 300 1 -1   -1   150 $wFFBFSig      #Minimum at 1.16
							      	    
bsub -q $q -o $logdir/logkeys_32.log FitPdfSampleKEYS 32 100 1 1000 300 1 0.45 1.5   45 $wFFBF         #Minimum at 1.22
bsub -q $q -o $logdir/logkeys_33.log FitPdfSampleKEYS 33 100 1 1000 300 1 0.4  1.4   50 $wFFBF         #Minimum at 1.36
bsub -q $q -o $logdir/logkeys_34.log FitPdfSampleKEYS 34 100 1 1000 300 1 -1   1.3   60 $wFFBF         #Minimum at 0.51
bsub -q $q -o $logdir/logkeys_35.log FitPdfSampleKEYS 35 100 1 1000 300 1 -1   1.5   45 $wFFBF         #Minimum at 0.58
bsub -q $q -o $logdir/logkeys_36.log FitPdfSampleKEYS 36  70 1 1000 300 1 -1   1.4   40 $wFFBF         #Minimum at 0.94

bsub -q $q -o $logdir/logkeys_41.log FitPdfSampleKEYS 41 110 1 1000 300 1 1    1.85  60 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_42.log FitPdfSampleKEYS 42 110 1 1000 300 1 0.8  1.95  45 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_43.log FitPdfSampleKEYS 43 115 1 1000 300 1 0.85 1.9   55 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_44.log FitPdfSampleKEYS 44 110 1 1000 300 1 0.75 1.9   55 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_45.log FitPdfSampleKEYS 45 130 1 1000 300 1 0.4  1.45  60 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_46.log FitPdfSampleKEYS 46 120 1 1000 300 1 0.45 1.5   55 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_47.log FitPdfSampleKEYS 47 120 1 1000 300 1 0.4  1.5   60 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_48.log FitPdfSampleKEYS 48 100 1 1000 300 1 0.4  1.6   55 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_51.log FitPdfSampleKEYS 51 150 1 1000 300 1 0.8  1.95  60 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_52.log FitPdfSampleKEYS 52 110 1 1000 300 1 0.8  1.95  60 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_53.log FitPdfSampleKEYS 53 130 1 1000 300 1 0.85 2.1   70 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_54.log FitPdfSampleKEYS 54 140 1 1000 300 1 0.75 2.05  70 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_55.log FitPdfSampleKEYS 55 130 1 1000 300 1 0.4  1.45  60 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_56.log FitPdfSampleKEYS 56 150 1 1000 300 1 0.45 1.5   60 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_57.log FitPdfSampleKEYS 57 160 1 1000 300 1 0.4  1.5   65 $wFFBF         #Minimum at -
bsub -q $q -o $logdir/logkeys_58.log FitPdfSampleKEYS 58 140 1 1000 300 1 0.4  1.6   60 $wFFBF         #Minimum at -
							     
bsub -q $q -o $logdir/logkeys_61.log FitPdfSampleKEYS 61  60 2.5 1000 300 1 0.8  1.5   25 $wFFBFSig      #Minimum at -
bsub -q $q -o $logdir/logkeys_62.log FitPdfSampleKEYS 62  50 2.5 1000 300 1 0.7  1.5   25 $wFFBFSig      #Minimum at -
bsub -q $q -o $logdir/logkeys_63.log FitPdfSampleKEYS 63  70 2.5 1000 300 1 0.7  1.5   30 $wFFBFSig      #Minimum at -
bsub -q $q -o $logdir/logkeys_64.log FitPdfSampleKEYS 64  65 2.5 1000 300 1 0.7  1.5   25 $wFFBFSig      #Minimum at -
bsub -q $q -o $logdir/logkeys_71.log FitPdfSampleKEYS 71  95 1 1000 300 1 0.8  1.6   50 $wFFBFSig      #Minimum at -
bsub -q $q -o $logdir/logkeys_72.log FitPdfSampleKEYS 72  85 1 1000 300 1 0.7  1.65  50 $wFFBFSig      #Minimum at -
bsub -q $q -o $logdir/logkeys_73.log FitPdfSampleKEYS 73  90 1 1000 300 1 0.55 1.5   55 $wFFBFSig      #Minimum at -
bsub -q $q -o $logdir/logkeys_74.log FitPdfSampleKEYS 74 120 1 1000 300 1 0.6  1.7   60 $wFFBFSig      #Minimum at -


bsub -q $q -o $logdir/logkeys_61.log FitPdfSampleKEYS 61  40 350 1000 300 30 0.4  1.9  220 $wFFBFSig      #Minimum at -, 61 
bsub -q $q -o $logdir/logkeys_62.log FitPdfSampleKEYS 62  40 350 1000 300 30 0.4  1.9  220 $wFFBFSig      #Minimum at -, 62
bsub -q $q -o $logdir/logkeys_63.log FitPdfSampleKEYS 63  40 350 1000 300 30 0.4  1.9  220 $wFFBFSig      #Minimum at -, 63
bsub -q $q -o $logdir/logkeys_64.log FitPdfSampleKEYS 64  40 350 1000 300 30 0.4  1.9  220 $wFFBFSig      #Minimum at -, 64
bsub -q $q -o $logdir/logkeys_71.log FitPdfSampleKEYS 71  55 350 1000 300 30 0.4  1.9  220 $wFFBFSig      #Minimum at -, 71
bsub -q $q -o $logdir/logkeys_72.log FitPdfSampleKEYS 72  55 350 1000 300 30 0.4  1.9  220 $wFFBFSig      #Minimum at -, 72
bsub -q $q -o $logdir/logkeys_73.log FitPdfSampleKEYS 73  55 350 1000 300 30 0.4  1.9  220 $wFFBFSig      #Minimum at -, 73
bsub -q $q -o $logdir/logkeys_74.log FitPdfSampleKEYS 74  55 350 1000 300 30 0.4  1.9  220 $wFFBFSig      #Minimum at -, 74


#bsub -q $q -o $logdir/logkeys_37.log FitPdfSampleKEYS 37 100 a 700 200 1 0.4  1.4   50 $wFFBF         #Minimum at 1.36

#set logdir = AWG82/logfiles/results/keys/Times2



#bsub -q $q -o $logdir/logkeysNoA${i}_70.log bbrroot -b -q 'babar_code/keysFit/keys.C("'${i}'","70",1,"NoR","nv",800,400)'
#while ( $n <= 31 )
#	rm $logdir/logkeys_${n}.log
#	set w = babar_code/Reweight/wAllFF.txt
#	if ($n == 6 || $n == 7 || $n == 16 || $n == 17 || $n >= 20) set w = wFFBFSig
#	bsub -q $q -o $logdir/logkeys_${n}.log FitPdfSampleKEYS $n 140 a 1000 300 1 -1 -1 $wFFBF
#	@ n = $n + 1
#end
