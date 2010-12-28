######################################################################
# Definitions
# $< is the suffix rule only variable of the trigger (file.c)
# $* is the suffix rule only variable of the prefix  (file)
######################################################################
FC=gfortran
CFLAG=-c
OFLAG=-o
INSTALL_DIR=/IPNShome/pfpeterson/bin
INSTALL_MOD=755

LIBFILES = dsort.o tred2.o tqli.o cgmin.o ran1.o

######################################################################
# target : dependencies command
######################################################################

all: allsig allppdf lendis convolute

allsig : allsig.o $(LIBFILES)
	$(FC) allsig.o $(LIBFILES) $(OFLAG) $@

allppdf : allppdf.f
	$(FC) $< $(OFLAG) $@

convolute : convolute.f
	$(FC) $< $(OFLAG) $@

lendis : lendis.f
	$(FC) $< $(OFLAG) $@

install: allsig allppdf lendis convolute
	cp allsig $(INSTALL_DIR)
	chmod $(INSTALL_MOD) $(INSTALL_DIR)/allsig
	cp allppdf $(INSTALL_DIR)
	chmod $(INSTALL_MOD) $(INSTALL_DIR)/allppdf
	cp lendis $(INSTALL_DIR)
	chmod $(INSTALL_MOD) $(INSTALL_DIR)/lendis
	cp convolute $(INSTALL_DIR)
	chmod $(INSTALL_MOD) $(INSTALL_DIR)/convolute


clean :	
	rm -f *.o allsig allppdf lendis convolute
