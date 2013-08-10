  proc donutGoodTracks { sequence DonutGoodTracksInput DonutGoodTracksOutput } {
     mod clone SmpMergerDefiner $DonutGoodTracksOutput
     seq append $sequence $DonutGoodTracksOutput
     catch { setProduction $DonutGoodTracksOutput }
     talkto $DonutGoodTracksOutput {
         inputListNames    set GoodTracksVeryLoose
         inputListNames    set V0BtaTrack
     }
  }

