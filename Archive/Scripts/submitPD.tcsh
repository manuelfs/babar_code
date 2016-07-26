#!/bin/tcsh

cd ..
srtpath `cat .current` $BFARCH
cd -

set q = medium
set d = `date +%y-%m-%d`
set logdir = AWG82/logfiles/results/keys/PDFitPdf/$d/
mkdir -p $logdir

set sample = RAll

FitPDpdf  1   130   30   1.65 
FitPDpdf  2   130   30   1.65 
FitPDpdf  3   130   30   1.65 
FitPDpdf  4   130   30   1.65 
FitPDpdf  5   130   30   1.65 
FitPDpdf  6   130   30   1.65 
FitPDpdf  7   130   30   1.65 
FitPDpdf  8   130   30   1.65 
FitPDpdf  9   130   30   1.65 
FitPDpdf 10   130   30   1.65 
FitPDpdf 11   130   30   1.65 
FitPDpdf 12   130   30   1.65 
FitPDpdf 13   130   30   1.65 
FitPDpdf 14   130   30   1.65 
FitPDpdf 15   130   30   1.65 
FitPDpdf 16   130   30   1.65 
FitPDpdf 17   130   30   1.65 
FitPDpdf 18   130   30   1.65 
FitPDpdf 19   130   30   1.65 
FitPDpdf 20   130   30   1.65 
FitPDpdf 41   130   30   1.65 
FitPDpdf 42   130   30   1.65 
FitPDpdf 43   130   30   1.65 
FitPDpdf 44   130   30   1.65 
FitPDpdf 51   130   30   1.65 
FitPDpdf 52   130   30   1.65 
FitPDpdf 53   130   30   1.65 
FitPDpdf 54   130   30   1.65 
FitPDpdf 61   130   30   1.65 
FitPDpdf 62   130   30   1.65 
FitPDpdf 63   130   30   1.65 
FitPDpdf 64   130   30   1.65 
					

echo "1.5   13" > ! babar_code/Systematics/MmissCut.txt
cpf keys/root/PDfitFinal_150_1300/* keys/root/PDfitFinal/
rut
.x babar_code/Systematics/YieldsPD.C+
.q
PDFit FitAll/PD/configData.txt 0 1


echo "1.5   5" > ! babar_code/Systematics/MmissCut.txt
cpf keys/root/PDfitFinal_150_500/* keys/root/PDfitFinal/
rut
.x babar_code/Systematics/YieldsPD.C+
.q
PDFit FitAll/PD/configData.txt 0 1

echo "5   13" > ! babar_code/Systematics/MmissCut.txt
cpf keys/root/PDfitFinal_500_1300/* keys/root/PDfitFinal/
rut
.x babar_code/Systematics/YieldsPD.C+
.q
PDFit FitAll/PD/configData.txt 0 1


#bsub -q $q -o $logdir/logFit_01.log FitPDpdf  1   130   30   1.65 # Sample 01
#bsub -q $q -o $logdir/logFit_02.log FitPDpdf  2   130   30   1.65 # Sample 02
#bsub -q $q -o $logdir/logFit_03.log FitPDpdf  3   130   30   1.65 # Sample 03
#bsub -q $q -o $logdir/logFit_04.log FitPDpdf  4   130   30   1.65 # Sample 04
#bsub -q $q -o $logdir/logFit_05.log FitPDpdf  5   130   30   1.65 # Sample 05
#bsub -q $q -o $logdir/logFit_06.log FitPDpdf  6   130   30   1.65 # Sample 06
#bsub -q $q -o $logdir/logFit_07.log FitPDpdf  7   130   30   1.65 # Sample 07
#bsub -q $q -o $logdir/logFit_08.log FitPDpdf  8   130   30   1.65 # Sample 08
#bsub -q $q -o $logdir/logFit_09.log FitPDpdf  9   130   30   1.65 # Sample 09
#bsub -q $q -o $logdir/logFit_10.log FitPDpdf 10   130   30   1.65 # Sample 10
#bsub -q $q -o $logdir/logFit_11.log FitPDpdf 11   130   30   1.65 # Sample 11
#bsub -q $q -o $logdir/logFit_12.log FitPDpdf 12   130   30   1.65 # Sample 12
#bsub -q $q -o $logdir/logFit_13.log FitPDpdf 13   130   30   1.65 # Sample 13
#bsub -q $q -o $logdir/logFit_14.log FitPDpdf 14   130   30   1.65 # Sample 14
#bsub -q $q -o $logdir/logFit_15.log FitPDpdf 15   130   30   1.65 # Sample 15
#bsub -q $q -o $logdir/logFit_16.log FitPDpdf 16   130   30   1.65 # Sample 16
#bsub -q $q -o $logdir/logFit_17.log FitPDpdf 17   130   30   1.65 # Sample 17
#bsub -q $q -o $logdir/logFit_18.log FitPDpdf 18   130   30   1.65 # Sample 18
#bsub -q $q -o $logdir/logFit_19.log FitPDpdf 19   130   30   1.65 # Sample 19
#bsub -q $q -o $logdir/logFit_20.log FitPDpdf 20   130   30   1.65 # Sample 20
#bsub -q $q -o $logdir/logFit_41.log FitPDpdf 41   130   30   1.65 # Sample 41
#bsub -q $q -o $logdir/logFit_42.log FitPDpdf 42   130   30   1.65 # Sample 42
#bsub -q $q -o $logdir/logFit_43.log FitPDpdf 43   130   30   1.65 # Sample 43
#bsub -q $q -o $logdir/logFit_44.log FitPDpdf 44   130   30   1.65 # Sample 44
#bsub -q $q -o $logdir/logFit_51.log FitPDpdf 51   130   30   1.65 # Sample 51
#bsub -q $q -o $logdir/logFit_52.log FitPDpdf 52   130   30   1.65 # Sample 52
#bsub -q $q -o $logdir/logFit_53.log FitPDpdf 53   130   30   1.65 # Sample 53
#bsub -q $q -o $logdir/logFit_54.log FitPDpdf 54   130   30   1.65 # Sample 54
#bsub -q $q -o $logdir/logFit_61.log FitPDpdf 61   130   30   1.65 # Sample 61
#bsub -q $q -o $logdir/logFit_62.log FitPDpdf 62   130   30   1.65 # Sample 62
#bsub -q $q -o $logdir/logFit_63.log FitPDpdf 63   130   30   1.65 # Sample 63
#bsub -q $q -o $logdir/logFit_64.log FitPDpdf 64   130   30   1.65 # Sample 64
							   
















































































































 


























































