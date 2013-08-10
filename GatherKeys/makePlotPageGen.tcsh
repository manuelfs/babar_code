#!/bin/tcsh

# Creates a Latex page with plots (from list of .eps files)
# argv[1] = output file name (for tex and ps. Tag .ps and .tex will be added)

set d = `date +%y-%m-%d`
if ( $# < 1 ) then
    set Tag = All
else
    set Tag = $1
endif

set File = ${d}_$Tag
set outfile_tex = ${File}.tex
set outfile_ps = ${File}.ps
set outfile_dvi = ${File}
set Plotfolder = /u/br/manuelf/releases/mike50/workdir/keys/eps/Plot/

rm $outfile_tex
cat template_begin.tex > $outfile_tex

# loop over list of plots
@ iplot = 0
foreach cplot ( "Fit  0: \(\\Dz\\tau\\nu\\Rightarrow\\Dz\)              " \
    "Fit  1: \(\\Dstarz\\tau\\nu\\Rightarrow\\Dstarz\)                  " \
    "Fit  2: \(\\Dp\\tau\\nu\\Rightarrow\\Dp\)                          " \
    "Fit  3: \(\\Dstarp\\tau\\nu\\Rightarrow\\Dstarp\)                  " \
    "Fit  8: \(\\Dstarz\\tau\\nu\\Rightarrow\\Dz\)                      " \
    "Fit  9: \(\\Dz\\tau\\nu\\Rightarrow\\Dstarz\)                      " \
    "Fit 10: \(\\Dstarp\\tau\\nu\\Rightarrow\\Dp\)                      " \
    "Fit 11: \(\\Dp\\tau\\nu\\Rightarrow\\Dstarp\)                      " \
    "Fit 12: \(\\Dz\\ell\\nu\\Rightarrow\\Dz\)                          " \
    "Fit 13: \(\\Dstarz\\ell\\nu\\Rightarrow\\Dstarz\)                  " \
    "Fit 14: \(\\Dp\\ell\\nu\\Rightarrow\\Dp\)                          " \
    "Fit 15: \(\\Dstarp\\ell\\nu\\Rightarrow\\Dstarp\)                  " \
    "Fit 20: \(\\Dstarz\\ell\\nu\\Rightarrow\\Dz\)                      " \
    "Fit 21: \(\\Dz\\ell\\nu\\Rightarrow\\Dstarz\)                      " \
    "Fit 22: \(\\Dstarp\\ell\\nu\\Rightarrow\\Dp\)                      " \
    "Fit 23: \(\\Dp\\ell\\nu\\Rightarrow\\Dstarp\)                      " \
    "Fit 26: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dz\)             " \
    "Fit 27: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\)         " \
    "Fit 28: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dp\)             " \
    "Fit 29: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\)         " \
    "Fit 30: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dz\\piz\)        " \
    "Fit 31: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\\piz\)    " \
    "Fit 32: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dp\\piz\)        " \
    "Fit 33: \(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\\piz\)    " \
    "Fit 38: \(\\Dstarz(\\ell/\\tau)\\nu\\Rightarrow\\Dz\\piz\)         " \
    "Fit 39: \(\\Dstarz(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\\piz\)     " \
    "Fit 40: \(\\Dstarp(\\ell/\\tau)\\nu\\Rightarrow\\Dp\\piz\)         " \
    "Fit 41: \(\\Dstarp(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\\piz\)     " \
    "Fit 42: \(\\Dz(\\ell/\\tau)\\nu\\Rightarrow\\Dz\\piz\)             " \
    "Fit 43: \(\\Dz(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\\piz\)         " \
    "Fit 44: \(\\Dp(\\ell/\\tau)\\nu\\Rightarrow\\Dp\\piz\)             " \
    "Fit 45: \(\\Dp(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\\piz\)         " \
    "Fit 46: Comb. background, signal                                   " \
    "Fit 47: Comb. background, \(\\Dpsp\\pi^0\)                         " \
    "Fit 48: Crossfeed in \(D\) channels, signal                        " \
    "Fit 49: Crossfeed in \(\\Dstar\) channels, signal                  " \
    "Fit 50: Crossfeed in \(\\Dpsp\\pi^0\)                              " \
    "Fit 51: Comb. background, \(\\Dz\)                                 " \
    "Fit 52: Comb. background, \(\\Dstarz\)                             " \
    "Fit 53: Comb. background, \(\\Dp\)                                 " \
    "Fit 54: Comb. background, \(\\Dstarp\)                             " \
    "Fit 55: Crossfeed in \(\\Dz\)                                      " \
    "Fit 56: Crossfeed in \(\\Dstarz\)                                  " \
    "Fit 57: Crossfeed in \(\\Dp\)                                      " \
    "Fit 58: Crossfeed in \(\\Dstarp\)                                  " \
    "Fit 59: Comb. background, \(\\Dz\\pi^0\)                           " \
    "Fit 60: Comb. background, \(\\Dstarz\\pi^0\)                       " \
    "Fit 61: Comb. background, \(\\Dp\\pi^0\)                           " \
    "Fit 62: Comb. background, \(\\Dstarp\\pi^0\)                       " \
    "Fit 63: Crossfeed in \(\\Dz\\pi^0\)                                " \
    "Fit 64: Crossfeed in \(\\Dstarz\\pi^0\)                            " \
    "Fit 65: Crossfeed in \(\\Dp\\pi^0\)                                " \
    "Fit 66: Crossfeed in \(\\Dstarp\\pi^0\)                            " )

  if ( $iplot == 4 ) @ iplot = 8
  if ( $iplot == 16 ) @ iplot = 20
  if ( $iplot == 24 ) @ iplot = 26
  if ( $iplot == 34 ) @ iplot = 38
  cat template_1x1.tex >> $outfile_tex
  set placeholder = "plot.eps"
  set plot = ${Plotfolder}${Tag}epsKeys_${iplot}_Fit.eps
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


