#!/usr/local/bin/tcsh

set wFFBFSig = babar_code/Reweight/wFFBFSig.txt
set wFFBF = babar_code/Reweight/wFFBF.txt
set q = xlong
set x = xxl


#echo Copying all fits to keys/root/Fit
#cpf keys/root/fitSep/* keys/root/Fit/

set d = `date +%y-%m-%d`
set sample = RestsnoSep
set logdir = AWG82/logfiles/results/keys/single/

mkdir -p $logdir
set logfile = $logdir/logGenFit_${sample}_${d}.log
rm $logfile


FitGenSample keys/root/fitSep/pdfKeys_12_Fit.root keys/root/fitSep/pdfKeys_16_Fit.root $sample keys/root/fitSep/pdfKeys_86_Fit.root > $logfile
FitGenSample keys/root/fitSep/pdfKeys_13_Fit.root keys/root/fitSep/pdfKeys_17_Fit.root $sample keys/root/fitSep/pdfKeys_86_Fit.root >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_14_Fit.root keys/root/fitSep/pdfKeys_18_Fit.root $sample keys/root/fitSep/pdfKeys_86_Fit.root >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_15_Fit.root keys/root/fitSep/pdfKeys_19_Fit.root $sample keys/root/fitSep/pdfKeys_86_Fit.root >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_20_Fit.root keys/root/fitSep/pdfKeys_24_Fit.root $sample keys/root/fitSep/pdfKeys_87_Fit.root >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_22_Fit.root keys/root/fitSep/pdfKeys_25_Fit.root $sample keys/root/fitSep/pdfKeys_87_Fit.root >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_26_Fit.root keys/root/fitSep/pdfKeys_70_Fit.root $sample >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_27_Fit.root keys/root/fitSep/pdfKeys_71_Fit.root $sample >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_28_Fit.root keys/root/fitSep/pdfKeys_72_Fit.root $sample >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_29_Fit.root keys/root/fitSep/pdfKeys_73_Fit.root $sample >> $logfile	   
FitGenSample keys/root/fitSep/pdfKeys_30_Fit.root keys/root/fitSep/pdfKeys_34_Fit.root $sample keys/root/fitSep/pdfKeys_74_Fit.root  >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_31_Fit.root keys/root/fitSep/pdfKeys_35_Fit.root $sample keys/root/fitSep/pdfKeys_75_Fit.root  >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_32_Fit.root keys/root/fitSep/pdfKeys_36_Fit.root $sample keys/root/fitSep/pdfKeys_76_Fit.root  >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_33_Fit.root keys/root/fitSep/pdfKeys_37_Fit.root $sample keys/root/fitSep/pdfKeys_77_Fit.root  >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_38_Fit.root keys/root/fitSep/pdfKeys_78_Fit.root $sample >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_39_Fit.root keys/root/fitSep/pdfKeys_79_Fit.root $sample >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_40_Fit.root keys/root/fitSep/pdfKeys_80_Fit.root $sample >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_41_Fit.root keys/root/fitSep/pdfKeys_81_Fit.root $sample >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_42_Fit.root keys/root/fitSep/pdfKeys_82_Fit.root $sample >> $logfile
FitGenSample keys/root/fitSep/pdfKeys_44_Fit.root keys/root/fitSep/pdfKeys_84_Fit.root $sample >> $logfile
											   

#FitGenSample keys/root/fitSep/pdfKeys_0_Fit.root  keys/root/fitSep/pdfKeys_4_Fit.root  Rest2noSep >  $logfile
#FitGenSample keys/root/fitSep/pdfKeys_1_Fit.root  keys/root/fitSep/pdfKeys_5_Fit.root  Rest2noSep >> $logfile
#FitGenSample keys/root/fitSep/pdfKeys_2_Fit.root  keys/root/fitSep/pdfKeys_6_Fit.root  Rest2noSep >> $logfile
#FitGenSample keys/root/fitSep/pdfKeys_3_Fit.root  keys/root/fitSep/pdfKeys_7_Fit.root  Rest2noSep >> $logfile
#FitGenSample keys/root/fitSep/pdfKeys_43_Fit.root keys/root/fitSep/pdfKeys_83_Fit.root Rest2noSep >> $logfile
#FitGenSample keys/root/fitSep/pdfKeys_45_Fit.root keys/root/fitSep/pdfKeys_85_Fit.root Rests2_1_2 >> $logfile
#echo Copying new fits to the final fit folder
#cpf keys/root/GenFit/* keys/root/Fit/












