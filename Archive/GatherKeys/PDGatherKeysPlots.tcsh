#!/bin/tcsh

# Creates a Latex page with plots (from list of .eps files)
# argv[1] = output file name (for tex and ps. Tag .ps and .tex will be added)

set d = `date +%y-%m-%d`
if ( $# < 1 ) then
    set Tag = RAll
else
    set Tag = $1
endif

set File = ${d}_$Tag
set outfile_tex = ${File}.tex
set outfile_ps = ${File}.ps
set outfile_dvi = ${File}
set Plotfolder = /u/br/manuelf/releases/mike50/workdir/keys/eps/fitPdf/

rm $outfile_tex
cat template_begin.tex > $outfile_tex

# loop over list of plots
@ iplot = 1
foreach cplot ( "Fit  1: \(\\Dz\\tau\\nu\\Rightarrow\\Dz\)              " \
    "Fit  2: \(\\Dstarz\\tau\\nu\\Rightarrow\\Dstarz\)                  " \
    "Fit  3: \(\\Dp\\tau\\nu\\Rightarrow\\Dp\)                          " \
    "Fit  4: \(\\Dstarp\\tau\\nu\\Rightarrow\\Dstarp\)                  " \
    "Fit  5: \(\\Dstarz\\tau\\nu\\Rightarrow\\Dz\)                      " \
    "Fit  6: \(\\Dz\\tau\\nu\\Rightarrow\\Dstarz\)                      " \
    "Fit  7: \(\\Dstarp\\tau\\nu\\Rightarrow\\Dp\)                      " \
    "Fit  8: \(\\Dp\\tau\\nu\\Rightarrow\\Dstarp\)                      " \
    "Fit  9: \(\\Dz\\ell\\nu\\Rightarrow\\Dz\)                          " \
    "Fit 10: \(\\Dstarz\\ell\\nu\\Rightarrow\\Dstarz\)                  " \
    "Fit 11: \(\\Dp\\ell\\nu\\Rightarrow\\Dp\)                          " \
    "Fit 12: \(\\Dstarp\\ell\\nu\\Rightarrow\\Dstarp\)                  " \
    "Fit 13: \(\\Dstarz\\ell\\nu\\Rightarrow\\Dz\)                      " \
    "Fit 14: \(\\Dz\\ell\\nu\\Rightarrow\\Dstarz\)                      " \
    "Fit 15: \(\\Dstarp\\ell\\nu\\Rightarrow\\Dp\)                      " \
    "Fit 16: \(\\Dp\\ell\\nu\\Rightarrow\\Dstarp\)                      " \
    "Fit 17: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dz\)             " \
    "Fit 18: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\)         " \
    "Fit 19: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dp\)             " \
    "Fit 20: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\)         " \
    "Fit 41: Crossfeed in \(\\Dz\)                                      " \
    "Fit 42: Crossfeed in \(\\Dstarz\)                                  " \
    "Fit 43: Crossfeed in \(\\Dp\)                                      " \
    "Fit 44: Crossfeed in \(\\Dstarp\)                                  " \
    "Fit 51: Comb. background, \(\\Dz\)                                 " \
    "Fit 52: Comb. background, \(\\Dstarz\)                             " \
    "Fit 53: Comb. background, \(\\Dp\)                                 " \
    "Fit 54: Comb. background, \(\\Dstarp\)                             " \
    "Fit 61: Continuum, \(\\Dz\)                                        " \
    "Fit 62: Continuum, \(\\Dstarz\)                                    " \
    "Fit 63: Continuum, \(\\Dp\)                                        " \
    "Fit 64: Continuum, \(\\Dstarp\)                                    " )

  if ( $iplot == 21 ) @ iplot = 41
  if ( $iplot == 45 ) @ iplot = 51
  if ( $iplot == 55 ) @ iplot = 61
  cat template_1x1.tex >> $outfile_tex
  set placeholder = "plot.eps"
  set plot = ${Plotfolder}/PDFit_${iplot}.eps
  #if( $iplot <= 9 ) set plot = ${Plotfolder}${Tag}epsKeys_0${iplot}_Fit.eps
  replace $placeholder $plot -- $outfile_tex 
  echo ${cplot}"}\n\\end{center}" >> $outfile_tex 
  @ iplot++
end

echo "\\end{document}" >> $outfile_tex

latex $outfile_tex
dvips $outfile_dvi
# ps2pdf $outfile_ps
rm *.aux
rm *.log
rm *.dvi

echo "\ngv $outfile_ps\n"
