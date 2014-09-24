CAD	= /home/ivan/oldstuff/my_espresso
CADROOT	= /home/ivan
SIS = /home/ivan/Sis
SHELL = /bin/sh
MAKE = /usr/bin/make
REQUIRE = 

DIRS = util sparse mincov pla

DIRD = 

TARGET	= my_espresso

TARLIB	= libmy_espresso.a

LIBSS 	= -lm -lpla -lsparse -lmincov -lutil
LIBSD	= 	  

CC	= gcc
CFLAGSS = -g -Wall
LDFLAGS = -g -L$(CAD)/lib 

#---------------------------------------------------------------------------

$(TARGET): $(TARLIB) do_link

$(TARLIB): do_make do_lns

$(TARLIBG): do_touch do_makeg do_lns
	mv $(TARLIB) $(TARLIBG)
$(TARLINT): do_lintlib do_repack_lint

do_link:
	$(CC) -o $(TARGET) main.c -I$(CAD)/include $(LDFLAGS) $(LIBSD) $(LIBSS) 

do_make:
	@for dir in $(DIRD); do			\
	    (cd $$dir; 					\
	     echo Making Dynamic $$dir ...;		\
	     ${MAKE} CC=$(CC) 'CFLAGS=$(CFLAGSD) $$(INCLUDE)' SIS=$(SIS) CADROOT=$(CADROOT) CAD=$(CAD) lib$$dir.so)\
	done
	@for dir in $(DIRS); do			\
	    (cd $$dir; 					\
	     echo Making Static $$dir ...;		\
	     ${MAKE} CC=$(CC) 'CFLAGS=$(CFLAGSS) $$(INCLUDE)' SIS=$(SIS) CADROOT=$(CADROOT) CAD=$(CAD)	lib$$dir.a)\
	done
	
do_lns:
	@for dir in $(DIRD); do	   \
	   (echo Making library link $$dir ...; \
	   ln -f -s $(CAD)/$$dir/lib$$dir.so $(CAD)/lib/lib$$dir.so) \
	done
	@for dir in $(DIRS); do	   \
	   (echo Making library link $$dir ...; \
	   ln -f -s $(CAD)/$$dir/lib$$dir.a $(CAD)/lib/lib$$dir.a) \
	done

clean: cleansome

cleansome:
	rm -f lib/*
	@for dir in $(DIRS); do				\
	    (cd $$dir; 					\
	     echo Cleaning $$dir ...; 			\
	     ${MAKE} -i CAD=$(CAD) SIS=$(SIS) strip_depend >/dev/null	\
	     ${MAKE} -i CAD=$(CAD) SIS=$(SIS) clean >/dev/null)	\
	done
	@for dir in $(DIRD); do				\
	    (cd $$dir; 					\
	     echo Cleaning $$dir ...; 			\
	     ${MAKE} -i CAD=$(CAD) SIS=$(SIS) strip_depend >/dev/null	\
	     ${MAKE} -i CAD=$(CAD) SIS=$(SIS) clean >/dev/null)	\
	done

backup:
	@rm -f my_espresso.tar
	@tar -cf my_espresso.tar Makefile main.c
	@for dir in $(DIRS); do			\
	  (tar -rf my_espresso.tar ./$$dir/*.c ./$$dir/*.h ./$$dir/Makefile)	\
        done					
	@for dir in $(DIRD); do 		\
          (tar -rf my_espresso.tar ./$$dir/*.c ./$$dir/*.h ./$$dir/Makefile)	\
	done
	@gzip -9 my_espresso.tar

	
wc:
	wc */*.[cly] | tail -1
	grep ';' */*.[cly] | wc | tail -1

compile:
	rm -f io/read_eqn.c io/eqnlex.c 
	rm -f genlib/readlib.c genlib/readliblex.c
	date; ${MAKE} $(TARGETG)
	date; ${MAKE} cleansome $(TARGET) 
	date; ${MAKE} $(TARLINT)
	date; ${MAKE} do_lint | grep -v '\.c:$$' | grep -v '^lint -I' >.lint.out 2>&1
	date; ${MAKE} tags
# hack -- adjust filenames for nfs access
#	sed 's!/users/awang!/net/beeblebrox/users/saldanha!g' tags >tags1
#	mv tags1 tags
#	${MAKE} depend
	${MAKE} help

require:
	@echo $(REQUIRE)

#-------------------------------------------------------------


# dummy targets
nothing:;

build header print uninstall debug debug-g debug-pg install.lint: nothing
	@echo no rule

_force:
