######################################################################
# Definitions
# $< is the suffix rule only variable of the trigger (file.c)
# $* is the suffix rule only variable of the prefix  (file)
######################################################################
F77   = f77
CFLAG = -c
OFLAG = -o
INSTALL_DIR = /IPNShome/pfpeterson/bin
INSTALL_MOD = 755

LIBFILES = dsort.o tred2.o tqli.o cgmin.o ran1.o

######################################################################
# dependency rools for suffixes
######################################################################
# to create an object file
.f.o : $(F77) $(CFLAG) $<
# to create an executable
.f: $(F77) $< $(OFLAG) $*

######################################################################
# target : dependencies command
######################################################################

all: allsig allppdf lendis convolute

allsig : allsig.o $(LIBFILES)
	$(F77) allsig.o $(LIBFILES) $(OFLAG) $@

allppdf : allppdf.f
	$(F77) $< $(OFLAG) $@

convolute : convolute.f
	$(F77) $< $(OFLAG) $@

lendis : lendis.f
	$(F77) $< $(OFLAG) $@

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
