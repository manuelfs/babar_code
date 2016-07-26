#!/usr/local/bin/tcsh

set wFFBFSig = babar_code/Reweight/wFFBFSig.txt
set wFFBF = babar_code/Reweight/wFFBF.txt
set q = long
set x = xlong


#echo Copying all fits to keys/root/Fit
#cpf keys/root/fitSep/* keys/root/Fit/

set d = `date +%y-%m-%d`
set sample = RestsnoSep
set logdir = AWG82/logfiles/results/keys/TimesGenFit/

mkdir -p $logdir
set logfile = $logdir/logGenFit_${sample}_${d}.log
rm $logfile
cpf keys/root/Times/hTimes_48_Fit.root keys/root/Times/hTimes_55_Fit.root 
cpf keys/root/Times/hTimes_48_Fit.root keys/root/Times/hTimes_57_Fit.root 
cpf keys/root/Times/hTimes_49_Fit.root keys/root/Times/hTimes_56_Fit.root 
cpf keys/root/Times/hTimes_49_Fit.root keys/root/Times/hTimes_58_Fit.root 
cpf keys/root/Times/hTimes_50_Fit.root keys/root/Times/hTimes_63_Fit.root 
cpf keys/root/Times/hTimes_50_Fit.root keys/root/Times/hTimes_64_Fit.root 
cpf keys/root/Times/hTimes_50_Fit.root keys/root/Times/hTimes_65_Fit.root 
cpf keys/root/Times/hTimes_50_Fit.root keys/root/Times/hTimes_66_Fit.root 

set dir = keys/root/Times
bsub -q $q -o $logdir/log_${sample}_${d}_12.log FitGenSample $dir/hTimes_12_Fit.root $dir/hTimes_16_Fit.root $sample 500 $dir/hTimes_86_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_13.log FitGenSample $dir/hTimes_13_Fit.root $dir/hTimes_17_Fit.root $sample 500 $dir/hTimes_86_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_14.log FitGenSample $dir/hTimes_14_Fit.root $dir/hTimes_18_Fit.root $sample 500 $dir/hTimes_86_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_15.log FitGenSample $dir/hTimes_15_Fit.root $dir/hTimes_19_Fit.root $sample 500 $dir/hTimes_86_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_20.log FitGenSample $dir/hTimes_20_Fit.root $dir/hTimes_24_Fit.root $sample 500 $dir/hTimes_87_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_22.log FitGenSample $dir/hTimes_22_Fit.root $dir/hTimes_25_Fit.root $sample 500 $dir/hTimes_87_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_26.log FitGenSample $dir/hTimes_26_Fit.root $dir/hTimes_70_Fit.root $sample 500 
bsub -q $q -o $logdir/log_${sample}_${d}_27.log FitGenSample $dir/hTimes_27_Fit.root $dir/hTimes_71_Fit.root $sample 500 
bsub -q $q -o $logdir/log_${sample}_${d}_28.log FitGenSample $dir/hTimes_28_Fit.root $dir/hTimes_72_Fit.root $sample 500 
bsub -q $q -o $logdir/log_${sample}_${d}_29.log FitGenSample $dir/hTimes_29_Fit.root $dir/hTimes_73_Fit.root $sample 500 	   
bsub -q $q -o $logdir/log_${sample}_${d}_30.log FitGenSample $dir/hTimes_30_Fit.root $dir/hTimes_34_Fit.root $sample 500 $dir/hTimes_74_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_31.log FitGenSample $dir/hTimes_31_Fit.root $dir/hTimes_35_Fit.root $sample 500 $dir/hTimes_75_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_32.log FitGenSample $dir/hTimes_32_Fit.root $dir/hTimes_36_Fit.root $sample 500 $dir/hTimes_76_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_33.log FitGenSample $dir/hTimes_33_Fit.root $dir/hTimes_37_Fit.root $sample 500 $dir/hTimes_77_Fit.root 
bsub -q $q -o $logdir/log_${sample}_${d}_38.log FitGenSample $dir/hTimes_38_Fit.root $dir/hTimes_78_Fit.root $sample 500 
bsub -q $q -o $logdir/log_${sample}_${d}_39.log FitGenSample $dir/hTimes_39_Fit.root $dir/hTimes_79_Fit.root $sample 500 
bsub -q $q -o $logdir/log_${sample}_${d}_40.log FitGenSample $dir/hTimes_40_Fit.root $dir/hTimes_80_Fit.root $sample 500 
bsub -q $q -o $logdir/log_${sample}_${d}_41.log FitGenSample $dir/hTimes_41_Fit.root $dir/hTimes_81_Fit.root $sample 500 
bsub -q $q -o $logdir/log_${sample}_${d}_42.log FitGenSample $dir/hTimes_42_Fit.root $dir/hTimes_82_Fit.root $sample 500 
bsub -q $q -o $logdir/log_${sample}_${d}_44.log FitGenSample $dir/hTimes_44_Fit.root $dir/hTimes_84_Fit.root $sample 500 
											   
mv keys/root/Times/hTimes_12_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_13_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_14_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_15_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_20_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_22_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_26_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_27_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_28_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_29_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_30_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_31_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_32_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_33_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_38_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_39_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_40_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_41_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_42_Fit.root keys/root/TimesTemp/
mv keys/root/Times/hTimes_44_Fit.root keys/root/TimesTemp/
											   


#FitGenSample keys/root/fitSep/pdfKeys_0_Fit.root  keys/root/fitSep/pdfKeys_4_Fit.root  Rest2noSep >  $logfile
#FitGenSample keys/root/fitSep/pdfKeys_1_Fit.root  keys/root/fitSep/pdfKeys_5_Fit.root  Rest2noSep >> $logfile
#FitGenSample keys/root/fitSep/pdfKeys_2_Fit.root  keys/root/fitSep/pdfKeys_6_Fit.root  Rest2noSep >> $logfile
#FitGenSample keys/root/fitSep/pdfKeys_3_Fit.root  keys/root/fitSep/pdfKeys_7_Fit.root  Rest2noSep >> $logfile
#FitGenSample keys/root/fitSep/pdfKeys_43_Fit.root keys/root/fitSep/pdfKeys_83_Fit.root Rest2noSep >> $logfile
#FitGenSample keys/root/fitSep/pdfKeys_45_Fit.root keys/root/fitSep/pdfKeys_85_Fit.root Rests2_1_2 >> $logfile
#echo Copying new fits to the final fit folder
#cpf keys/root/GenFit/* keys/root/Fit/












