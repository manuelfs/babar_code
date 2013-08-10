#!/usr/local/bin/tcsh

    cd ..
    srtpath `cat .current` $BFARCH
    anal50boot
    cd -
    printf "\n\n"

    
    CompareGen D            Rest2
    CompareGen Ds           Rest2
    CompareGen Dfu          Rest2
    CompareGen Dsfd         Rest2
    CompareGen Dtau         Rest2
    CompareGen Dstau        Rest2
    CompareGen Dstaufd      Rest2
    CompareGen Dssfdd       Rest2
    CompareGen Dssfd        Rest2
    CompareGen pi0DssD      Rest2
    CompareGen pi0DssDs     Rest2
    CompareGen pi0Dfu       Rest2
    CompareGen pi0Dsfu      Rest2

    CompareGen D            Rest
    CompareGen Ds           Rest
    CompareGen Dfu          Rest
    CompareGen Dsfd         Rest
    CompareGen Dtau         Rest
    CompareGen Dstau        Rest
    CompareGen Dstaufd      Rest
    CompareGen Dssfdd       Rest
    CompareGen Dssfd        Rest
    CompareGen pi0DssD      Rest
    CompareGen pi0DssDs     Rest
    CompareGen pi0Dfu       Rest
    CompareGen pi0Dsfu      Rest


