#!/usr/local/bin/tcsh

if ($#argv > 0) then
    set folder = $argv[1]
else
    set folder = keys/Archive/10-03-08/Fit/
endif

cpf ${folder}/* keys/root/Fit/

cpf keys/root/fitSep/pdfKeys_51_Fit.root keys/root/Fit/
cpf keys/root/fitSep/pdfKeys_52_Fit.root keys/root/Fit/
cpf keys/root/fitSep/pdfKeys_53_Fit.root keys/root/Fit/
cpf keys/root/fitSep/pdfKeys_54_Fit.root keys/root/Fit/
cpf keys/root/fitSep/pdfKeys_59_Fit.root keys/root/Fit/
cpf keys/root/fitSep/pdfKeys_60_Fit.root keys/root/Fit/
cpf keys/root/fitSep/pdfKeys_61_Fit.root keys/root/Fit/
cpf keys/root/fitSep/pdfKeys_62_Fit.root keys/root/Fit/

cpf ${folder}/pdfKeys_48_Fit.root keys/root/Fit/pdfKeys_55_Fit.root
cpf ${folder}/pdfKeys_48_Fit.root keys/root/Fit/pdfKeys_57_Fit.root
cpf ${folder}/pdfKeys_49_Fit.root keys/root/Fit/pdfKeys_56_Fit.root
cpf ${folder}/pdfKeys_49_Fit.root keys/root/Fit/pdfKeys_58_Fit.root
cpf ${folder}/pdfKeys_50_Fit.root keys/root/Fit/pdfKeys_63_Fit.root
cpf ${folder}/pdfKeys_50_Fit.root keys/root/Fit/pdfKeys_64_Fit.root
cpf ${folder}/pdfKeys_50_Fit.root keys/root/Fit/pdfKeys_65_Fit.root
cpf ${folder}/pdfKeys_50_Fit.root keys/root/Fit/pdfKeys_66_Fit.root

