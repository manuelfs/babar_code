proc donutGoodPhotons { sequence DonutGoodPhotonsInput DonutGoodPhotonsOutput } {
    
    module clone GoodPhotonsVisibleESelection DonutGoodPhotons${sequence}
    talkto DonutGoodPhotons${sequence} {
	inputBtaCandidateList  set $DonutGoodPhotonsInput
	outputBtaCandidateList set $DonutGoodPhotonsOutput
	outputMap              set $DonutGoodPhotonsOutput
	eRawCut                set 0.050
	nCrys                  set 3
	latShapeCut            set 0.6
	thetaMin               set 0.32
	thetaMax               set 2.44
	TrkListForDeltaAlpha   set ChargedTracks 
	electronSelector       set PidLHElectrons
	maxNearestTrkDistance  set 999
	deltaAlphaMin          set 0.08
	# activate histogramming for this module
	histoFlag              set f
	# specify the number of entries in the ntuple (-1 means no limit)
	maxEventTuple          set 0
	#activate debugging histograms
	debugHistograms        set f
    }
    sequence append $sequence DonutGoodPhotons$sequence
}
