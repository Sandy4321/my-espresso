#/*
# * Revision Control Information
# *
# * $Source: /vol/opua/opua2/sis/sis-1.2/common/src/sis/util/RCS/Makefile,v $
# * $Author: sis $
# * $Revision: 1.3 $
# * $Date: 22/.0/.1 .0:.1:.2 $
# *
# */
# util -- miscellaneous utility routines
#---------------------------------------------------------------------------

CAD	= /projects/octtools/octtools/$(MACHINE)
SIS	= ..
LINTCREATEFLAG = -C

CC	= gcc

P	= util
PSRC	= cpu_time.c cpu_stats.c getopt.c safe_mem.c strsav.c texpand.c \
	  prtime.c pipefork.c saveimage.c pathsearch.c stub.c \
	  restart.c state.c tmpfile.c
POBJ	= $(PSRC:.c=.o)
PHDR	= util.h
MISC	= $(P).doc Makefile test-restart.c

TARGET	= libutil.a
HDR	=

LIBS	=
LINTLIBS=
INCLUDE	= -I/home/i/my_espresso/include

FLAGS	= -DUNIX
CFLAGS	= -g $(INCLUDE) -Wall
LDFLAGS	= -g

.c.o:
	$(CC) -c $(CFLAGS) $(FLAGS) $*.c
#---------------------------------------------------------------------------

$(TARGET): $(PHDR) $(OBJ) $(POBJ) $(LIBS)
#	$(CC) $(LDFLAGS) -o $(TARGET) $(OBJ) $(POBJ) $(LIBS)


install: lib$(P).a llib-l$(P).ln

lib$(P).a: $(POBJ)
	ar cr lib$(P).a $(POBJ)
	ranlib lib$(P).a

unpack: lib$(P).so
	@for i in $(POBJ); do				\
	    ln -s $(SIS)/$(P)/$$i $(SIS)/unpack;	\
	done

llib-l$(P).ln: $(PSRC) $(PHDR)
	lint $(LINTFLAGS) ${LINTCREATEFLAG}$(P) -n $(PSRC)

clean:
	rm -f $(TARGET) *.so *.ln *.o \
	    [Tt]ags [Mm]ake.out lint malloc.out gmon.out __.SYMDEF

tags: _force
	@for i in $(PSRC) $(PHDR); do			\
	    cwd=`pwd`; ctags -a $$cwd/$$i;		\
	done;

strip_depend:
	sed '/^#--DO NOT CHANGE ANYTHING AFTER THIS LINE/,$$d' Makefile >mktemp
	mv mktemp Makefile

depend:
	sed '/^#--DO NOT CHANGE ANYTHING AFTER THIS LINE/,$$d' Makefile >mktemp
	echo '#--DO NOT CHANGE ANYTHING AFTER THIS LINE' >>mktemp
	$(CAD)/bin/cc-M $(CFLAGS) $(PSRC) | \
	    sed 's|$(CAD)|$$(CAD)|g' | \
	    grep -v '/usr/include' >>mktemp
	mv mktemp Makefile

#-------------------------- IBM 3090 support -----------------
IBMHOST	= beeblebrox
IBMDIST	= 
ibmdist: $(PSRC) $(PHDR)
	rdist -Richw $(PSRC) $(PHDR) $(IBMHOST):$(IBMDIST)
#-------------------------------------------------------------
_force:

