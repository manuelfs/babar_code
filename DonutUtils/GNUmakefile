# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := 
LIBTMPLFILES :=
LIBDDLORDERED :=


#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS :=

$(addprefix $(bindir),$(EXTRABINS)): $(bindir)% : %.o

#################################################################
#++ binary products				[build it with 'bin']

#	ConstraintWiggler NormPlots PRDNormPlots GenCombWeights EffTable \
#	Q2BlindApp Q2CS MesFits KsBlindApp KsCS GenBFWeights GenFFWeights \
#	DstarstarBlindApp CalcPenalties GenPi0Weights MesTest SemilepCS YieldsApp \
#CXXFLAGS += -pg 
#LDFLAGS += -pg 
BINS := Q2_Systematics PRD_Fit prdFit prdAfterFit PRL_Plot EventMixer MultiEventMixer GenPdfSamples DonutFitFinal  plotVariable PlotFit Systematics KeysPlot CSample \
	SidebandmES WeightVariable ResolutionDss ResolutionFit FitPdfSampleKEYS \
	WeightPl VaryConstraints CrossValidation OptiBDT BootstrapBands IterBF
BINCCFILES :=  $(BINS:=.cc) 
# $(EXTRABINS:=.cc)

#++ Binary rules		 [in production, build it with 'bin']

$(addprefix $(bindir),$(BINS)): $(bindir)% : %.o

#++ shell script products.. 			[build it with 'bin']
#BINSCRIPTS := mkskel

#################################################################
#++ regression test scripts			[build it with 'test']

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
