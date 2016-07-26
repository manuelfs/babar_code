#!/bin/tcsh

if ( $# < 2 ) then
    echo "The format is: sendChi2.tcsh mBins pBins [samples=all]"
    printf "\n\n"
    exit -1
endif
set mBins = $1
set pBins = $2
if ( $# > 2 ) then
    set AllSamples = $3
else
    set AllSamples = all
endif
if ( $AllSamples == all ) set AllSamples = (030 050 075 100 120 150 200 250 300)
foreach sample ( $AllSamples )
    PlotFit FitAll/fits/TextFinalDataNe2x${sample}.txt sig${sample} All candM2 $mBins 1.5 10 histoPl 1 1.5 12 0.85 $pBins
end
