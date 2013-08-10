# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile,v 1.6 2012/08/23 02:12:07 manuelf Exp $
#
# 20090505  Use $(wildcard...) to expand "~manuelf/..." path.  Use SYSLIB
#	    instead of LOADLIBES to position TMVA library correctly.
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := 
LIBTMPLFILES := 
LIBDDLORDERED := 

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS := 
#$(addprefix $(bindir),$(EXTRABINS)): $(bindir)% : %.o

#################################################################
#++ binary products				[build it with 'bin']

BINS := DonutUserApp

BINCCFILES := DonutUserAUB.cc

#++ Binary rules		 [in production, build it with 'bin']

$(bindir)DonutUserApp : DonutUserAUB.o

#++ shell script products.. 			[build it with 'bin']
BINSCRIPTS := 

#################################################################
#++ regression test scripts			[build it with 'test']


#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
