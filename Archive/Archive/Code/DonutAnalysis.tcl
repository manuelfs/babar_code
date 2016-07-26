#------------------------------------------------------------------------------
# $Id: DonutAnalysis.tcl,v 1.4 2009/09/03 01:16:21 manuelf Exp $
# Sample MyMiniAnalysis.tcl file
#------------------------------------------------------------------------------
# always source the error logger early in your main tcl script
sourceFoundFile ErrLogger/ErrLog.tcl
sourceFoundFile FrameScripts/FwkCfgVar.tcl
sourceFoundFile FrameScripts/talkto.tcl
# Disable the use of envvars
set ProdTclOnly true

# set the error logging level to 'warning'.  If you encounter a configuration 
# error you can get more information using 'trace'
ErrLoggingLevel warning

## allowed values of BetaMiniReadPersistence are (currently) "Kan", "Bdb"
##
FwkCfgVar BetaMiniReadPersistence Kan

## allowed (non-expert) values of levelOfDetail are "micro", "cache", "extend" 
## or "refit"
##
FwkCfgVar levelOfDetail "cache"

## allowed values of ConfigPatch are "Run1", "Run2" or "MC".  This MUST be set 
## consistent ## with your input data type or you will get INCONSISTENT OR 
## INCORRECT RESULTS
##
FwkCfgVar ConfigPatch "MC"
if [ regexp SP- $inputTcl ] {
    echo "It is Monte Carlo. Setting isData to true and ConfigPatch to MC"
    set isData false
    set ConfigPatch MC
} elseif [ regexp Run1 $inputTcl ] {
    echo "It is run1 data. Setting isData to false and ConfigPatch to Run1"
    set isData true
    set ConfigPatch Run1
} else {
    echo "It is run2-6 data. Setting isData to false and ConfigPatch to Run2"
    set isData true
    set ConfigPatch Run2
} 
#ErrMsg fatal "It's just a test. Exiting..."

##
## Set the number of events to run. If this isn't set, all events in the
## input collections will be processed.
##
FwkCfgVar NEvent 40

## choose the flavor of ntuple to write (hbook or root) and the file name
##
if [info exists rootFile] {
    set histFileName $rootFile
}
FwkCfgVar BetaMiniTuple "root"
FwkCfgVar histFileName "MyMiniAnalysis.root"

##
##  You can enter input collections two ways: either append them to a list, or
##  explicitly enter them in the input module. Do one or the other, BUT NOT 
##  BOTH.
##  If inputList is set before executing btaMini.tcl, that will automatically
##  add the collections to the appropriate input module, otherwise make sure you
##  talk to the right one.
##
## lappend inputList collection1 collection2 ...
##
##  OR THE FOLLOWING (choose the correct one based on persistence)
##
## talkto BdbEventInput {
## talkto KanEventInput {
##    input add collection1
##    input add collection2
##    ...
## }

if [info exists inputTcl] {
    sourceFoundFile $inputTcl
} else {
  ErrMsg fatal "Envvar inputTcl not set!"
}

sourceFoundFile BetaMiniUser/btaMini.tcl
#sourceFoundFile BetaMiniUser/btaMiniPhysics.tcl
#sourceFoundFile BetaMiniUser/btaMiniPhysProdSequence.tcl
sourceFoundFile UsrTools/UsrDataProcs.tcl
enableReadUsrData
readCandUsrData VcbCandData
sourceFoundFile PidTools/PidSequence.tcl
sequence append BetaMiniSequence PidSequence

sourceFoundFile SimpleComposition/SmpProcs.tcl
sourceFoundFile SimpleComposition/SmpPi0ProdSequence.tcl
sequence append BetaMiniSequence SmpPi0ProdSequence
sourceFoundFile CompositionSequences/CompPi0Sequence.tcl
sequence append BetaMiniSequence CompPi0Sequence
sequence enable CompPi0Sequence
# sourceFoundFile SimpleComposition/SmpComboProdSequence.tcl
# sequence append BetaMiniSequence SmpComboProdSequence
sourceFoundFile SimpleComposition/SmpKsProdSequence.tcl
sequence append BetaMiniSequence SmpKsProdSequence

if { $isData == "true" } {

} else {
    echo "manuelf: It's MC, setting PID tweaking and setting pepEnergiesCorrFile"
    pidCfg_dataset ascii
    set pidAsciiRevision r24c
    pidCfg_mode tweak BDTVeryLooseMuonSelection
    pidCfg_mode tweak NNVeryLooseMuonSelection
    pidCfg_mode tweak LooseKMElectronMicroSelection
    pidCfg_mode tweak PidLHElectronSelector
    pidCfg_mode tweak TightLHKaonMicroSelection
    pidCfg_mode tweak LooseLHPionMicroSelection

#    talkto PepBuildEnv {
#	set pepEnergiesCorrFile PepCond/pepEnergiesCorr.raw
#    }
}

sourceFoundFile DonutUser/DonutGoodTracks.tcl
sourceFoundFile DonutUser/DonutGoodPhotons.tcl

sequence create DonutSequence

donutGoodTracks  DonutSequence ChargedTracks RecoilGoodTracks
donutGoodPhotons DonutSequence CalorNeutral  RecoilGoodPhotons

talkto DonutGoodPhotonsDonutSequence {
    histoFlag set t
    debugHistograms set t
}
# talkto DonutGoodTracksDonutSequence {
#     histoFlag set t
#     debugHistograms set t
# }

sequence append DonutSequence DonutEVeto
sequence append DonutSequence DonutMCTagger
sequence append DonutSequence DonutAnalysis

talkto DonutMCTagger {
    Verbose set $isVerbose
    #Verbose set t
}
talkto DonutAnalysis {
    TrackList set RecoilGoodTracks
    PhotonList set RecoilGoodPhotons
    muonList set $MuonList
    elecList set $ElecList
    EmptyFile set $EmptymodesFile
    doBestBEEx set $BestBEEx
    doVertex set t
    doMass set f
    Verbose set $isVerbose
    deltaEcut set $dEcut
    mEScut set $mescut
    StemFolder set $Dfolder
    if [info exists env(SIDEBAND)] {
	echo "........................";
	echo "RUNNING IN SIDEBAND MODE";
	echo "........................";
	doSideband set t
    }
    if [info exists env(DELTAE)] {
	echo "........................";
	echo "RUNNING IN DELTAE SIDEBAND MODE";
	echo "........................";
	deltaESideband set t
    }
    if [info exists env(SMEAR)] {
	echo "........................";
	echo "RUNNING IN SMEAR MODE";
	echo "........................";
	doSmear set t
    }
    if [info exists env(SKIPBPLUS)] {
	echo "........................";
	echo "RUNNING IN SKIPBPLUS MODE";
	echo "........................";
	skipBPlus set t
    }
}

path append Everything DonutSequence

action enable SumTimeAction
#path list
echo "manuelf: pidAsciiRevision " $pidAsciiRevision
talkto EvtCounter {
    printFreq set $FreqPrint
#  skipStart set 1
#  skipStop  set 2040
}

if [info exists env(NORUN)] {

} else {
 if { $NEvent == "All" } {
     ev begin
     exit
 } elseif { $NEvent == "NoRun" } {
     echo "Returning control"
 } else {
     ev begin -nev $NEvent
     exit
 }
}
