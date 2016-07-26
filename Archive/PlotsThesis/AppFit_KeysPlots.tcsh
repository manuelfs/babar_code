#!/bin/tcsh

# Creates a Latex page with plots (from list of .eps files)

set File = public_html/AppFit_KeysPlots
set outfile_tex = ${File}.tex
set Plotfolder = images/KeysPlots/

rm $outfile_tex

# loop over list of plots
@ iplot = 1
foreach cplot ( "\(\\Dz\\tau\\nu\\Rightarrow\\Dz\)              " \
    "\(\\Dstarz\\tau\\nu\\Rightarrow\\Dstarz\)                  " \
    "\(\\Dp\\tau\\nu\\Rightarrow\\Dp\)                          " \
    "\(\\Dstarp\\tau\\nu\\Rightarrow\\Dstarp\)                  " \
    "\(\\Dstarz\\tau\\nu\\Rightarrow\\Dz\)                      " \
    "\(\\Dz\\tau\\nu\\Rightarrow\\Dstarz\)                      " \
    "\(\\Dstarp\\tau\\nu\\Rightarrow\\Dp\)                      " \
    "\(\\Dp\\tau\\nu\\Rightarrow\\Dstarp\)                      " \
    "\(\\Dz\\ell\\nu\\Rightarrow\\Dz\)                          " \
    "\(\\Dstarz\\ell\\nu\\Rightarrow\\Dstarz\)                  " \
    "\(\\Dp\\ell\\nu\\Rightarrow\\Dp\)                          " \
    "\(\\Dstarp\\ell\\nu\\Rightarrow\\Dstarp\)                  " \
    "\(\\Dstarz\\ell\\nu\\Rightarrow\\Dz\)                      " \
    "\(\\Dz\\ell\\nu\\Rightarrow\\Dstarz\)                      " \
    "\(\\Dstarp\\ell\\nu\\Rightarrow\\Dp\)                      " \
    "\(\\Dp\\ell\\nu\\Rightarrow\\Dstarp\)                      " \
    "\(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dz\)             " \
    "\(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\)         " \
    "\(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dp\)             " \
    "\(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\)         " \
    "\(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dz\\piz\)        " \
    "\(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\\piz\)    " \
    "\(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dp\\piz\)        " \
    "\(\\dstrstr(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\\piz\)    " \
    "\(\\Dstarz(\\ell/\\tau)\\nu\\Rightarrow\\Dz\\piz\)         " \
    "\(\\Dstarz(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\\piz\)     " \
    "\(\\Dstarp(\\ell/\\tau)\\nu\\Rightarrow\\Dp\\piz\)         " \
    "\(\\Dstarp(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\\piz\)     " \
    "\(\\Dz(\\ell/\\tau)\\nu\\Rightarrow\\Dz\\piz\)         " \
    "\(\\Dz(\\ell/\\tau)\\nu\\Rightarrow\\Dstarz\\piz\)     " \
    "\(\\Dp(\\ell/\\tau)\\nu\\Rightarrow\\Dp\\piz\)         " \
    "\(\\Dp(\\ell/\\tau)\\nu\\Rightarrow\\Dstarp\\piz\)     " \
    "Crossfeed in \(\\Dz\)                                      " \
    "Crossfeed in \(\\Dstarz\)                                  " \
    "Crossfeed in \(\\Dp\)                                      " \
    "Crossfeed in \(\\Dstarp\)                                  " \
    "Crossfeed in \(\\Dz\\pi^0\)                                " \
    "Crossfeed in \(\\Dstarz\\pi^0\)                            " \
    "Crossfeed in \(\\Dp\\pi^0\)                                " \
    "Crossfeed in \(\\Dstarp\\pi^0\)                            " \
    "Combinatorial background, \(\\Dz\)                                 " \
    "Combinatorial background, \(\\Dstarz\)                             " \
    "Combinatorial background, \(\\Dp\)                                 " \
    "Combinatorial background, \(\\Dstarp\)                             " \
    "Combinatorial background, \(\\Dz\\pi^0\)                           " \
    "Combinatorial background, \(\\Dstarz\\pi^0\)                       " \
    "Combinatorial background, \(\\Dp\\pi^0\)                           " \
    "Combinatorial background, \(\\Dstarp\\pi^0\)                       " \
    "Continuum, \(\\Dz\)                                        " \
    "Continuum, \(\\Dstarz\)                                    " \
    "Continuum, \(\\Dp\)                                        " \
    "Continuum, \(\\Dstarp\)                                    " \
    "Continuum, \(\\Dz\\pi^0\)                                  " \
    "Continuum, \(\\Dstarz\\pi^0\)                              " \
    "Continuum, \(\\Dp\\pi^0\)                                  " \
    "Continuum, \(\\Dstarp\\pi^0\)                              " )

  if ( $iplot == 33 ) @ iplot = 41
  if ( $iplot == 49 ) @ iplot = 51
  if ( $iplot == 59 ) @ iplot = 61
  cat babar_code/PlotsThesis/AppFit_KeysPlots_template.tex >> $outfile_tex
  set placeholder = "plot.eps"
  set plot = ${Plotfolder}RAllepsKeys_${iplot}_Fit.eps
  /opt/TWWfsw/bin/replace $placeholder $plot -- $outfile_tex 
  echo ${cplot}" component.}\n\\label{fig:AppFit_PDF"${iplot}"}\n\\end{figure}\n" >> $outfile_tex 
  @ iplot++
end

