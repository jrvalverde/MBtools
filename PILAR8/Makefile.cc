#   MAKEFILE for the Protein Identification and Library
#   Access Resource (PILAR)
#
#   Designed by:
#	J. R. Valverde	    (JRamon@emduam51.BitNet)
#			    (JRamon@ccuam1.uam.es)
#		      	    (JRamon@uamed.uam.es)

.SUFFIXES: .c .o

CC = cc

#   = -g3 -O3
#   = -c -g2 
#   = -c -g1 -O
#   = -DMEM_DEBUG	Add this flag if debugging memory allocation
# -Dconst=' ' -DUNIX -D__DATE__='"Jan 19, 1992"'
#   = -p -g -Ih -DMEM_DEBUG -D__STDC__ -Dconst=' ' -DUNIX -D__DATE__='"Jan 19, 1992"'
CFLAGS = -g -Ih -DMEM_DEBUG -D__STDC__ -Dconst=' ' -DUNIX -D__DATE__='"Sept 15, 1992"'


LDFLAGS	= -Llib -lPILAR -lm -lcurses -ltermcap

H	= h/portable.h h/P_defs.h h/P_types.h

PILIB	= lib/libPILAR.a

LIBS	= lib/libPlib.a lib/libPgen.a lib/libPdb.a lib/libPrg.a

# compile library modules, then strip unnecessary symbols
.c.o:
	$(CC) -c $(CFLAGS) $*.c -o $*.o
#	ld -r -x $*.o
#	mv a.out $*.o

#
#   MAKE PROGRAMS:
#	1. All: Make all the component programs
#

All: Pilar_2 Pilar_4 Pilar_5 Pilar_6 Pilar_7 Pilar_8 Pilar_9 Pilar_10

#
#	2. Dependencies for each program
#

Pilar_10: Pilar_10.c $H $(PILIB)
	cc $(CFLAGS) Pilar_10.c $(LDFLAGS) -o Pilar_10

Pilar_9: Pilar_9.c $H $(PILIB)
	cc $(CFLAGS) Pilar_9.c $(LDFLAGS) -o Pilar_9

Pilar_8: Pilar_8.c $H $(PILIB)
	cc $(CFLAGS) Pilar_8.c $(LDFLAGS) -o Pilar_8

Pilar_7: Pilar_7.c $H $(PILIB)
	cc $(CFLAGS) Pilar_7.c $(LDFLAGS) -o Pilar_7

Pilar_6: Pilar_6.c $H $(PILIB)
	cc $(CFLAGS) Pilar_6.c $(LDFLAGS) -o Pilar_6

Pilar_5: Pilar_5.c $H $(PILIB)
	cc $(CFLAGS) Pilar_5.c $(LDFLAGS) -o Pilar_5

Pilar_4: Pilar_4.c $H $(PILIB)
	cc $(CFLAGS) Pilar_4.c $(LDFLAGS) -o Pilar_4

Pilar_2: Pilar_2.c $H $(PILIB)
	cc $(CFLAGS) Pilar_2.c $(LDFLAGS) -o Pilar_2


#
#   MAKE LIBRARIES
#
#	1. PILAR specific routines
#

P_LIBS: P_lib P_gen P_db P_rg

SPECIFICS = mshell.o P_utils.o P_protein.o P_descr.o P_init.o P_region.o \
P_parser.o P_values.o P_think.o

P_lib: lib/libPlib.a

lib/libPlib.a: $(SPECIFICS)
	ar r lib/libPlib.a $?
	ranlib lib/libPlib.a

#   General utility routines

mshell.o: mshell.c h/mshell.h

P_utils.o: P_utils.c

#   Protein manipulation routines

P_protein.o: P_protein.c h/P_list.h h/P_db.h h/P_utils.h $H

P_init.o: P_init.c h/P_list.h $H

#   Structure prediction and manipulation routines

P_region.o: P_region.c h/P_region.h h/P_utils.h $H

P_values.o: P_values.c h/P_utils.h h/P_list.h $H

#   Rules manipulation routines

P_descr.o: P_descr.c h/P_parser.h h/P_list.h h/P_utils.h $H

P_parser.o: P_parser.c h/P_list.h h/P_utils.h $H

P_think.o: P_think.c h/P_list.h $H

#
#	2. Modules for manipulation of generic data structures:
#

GENERICS = P_AVL.o P_dlist.o P_gtree.o P_list.o P_queue.o \
P_sort.o P_stack.o P_stree.o

P_gen: lib/libPgen.a

lib/libPgen.a: $(GENERICS)
	ar r lib/libPgen.a $?
	ranlib lib/libPgen.a

P_AVL.o: P_AVL.c

P_dlist.o: P_dlist.c

P_gtree.o: P_gtree.c

P_list.o: P_list.c

P_queue.o: P_queue.c

P_sort.o: P_sort.c

P_stack.o: P_stack.c

P_stree.o: P_stree.c

#
#	3. Database access module
#
#	The DataBase Module is actually an object library made up of
#   one master submodule, P_DBANK, and several lower level modules that
#   perform the specific database access, one per supported data base
#   format.
#

DBFORMATS = db/P_EMBL.o db/P_Fitch.o db/P_GenBank.o db/P_GCGSw.o \
	db/P_IG.o db/P_Pearson.o db/P_PIR.o db/P_Preer.o db/P_Strider.o \
	db/P_unknown.o db/P_UWGCG.o db/P_Zuker.o

P_db: lib/libPdb.a

lib/libPdb.a: P_db.o $(DBFORMATS)
	ar r lib/libPdb.a $?
	ranlib lib/libPdb.a

P_db.o: P_db.c h/P_utils.h $H

db/P_EMBL.o: db/P_EMBL.c h/P_dblow.h $H

db/P_Fitch.o: db/P_Fitch.c h/P_dblow.h $H

db/P_GCGSw.o: db/P_GCGSw.c h/P_dblow.h $H

db/P_GenBank.o: db/P_GenBank.c h/P_dblow.h $H

db/P_IG.o: db/P_IG.c h/P_dblow.h $H

db/P_Pearson.o: db/P_Pearson.c h/P_dblow.h

db/P_PIR.o: db/P_PIR.c h/P_dblow.h $H

db/P_Preer.o: db/P_Preer.c h/P_dblow.h $H

db/P_Strider.o: db/P_Strider.c h/P_dblow.h $H

db/P_unknown.o: db/P_unknown.c h/P_dblow.h $H

db/P_UWGCG.o: db/P_UWGCG.c h/P_dblow.h $H

db/P_Zuker.o: db/P_Zuker.c h/P_dblow.h $H

#
#	4. Graphic libraries for ReGis mode
#

GRAPHICS = rg/P_graph1.o rg/P_graph2.o rg/P_graph3.o rg/P_graph4.o

P_rg: lib/libPrg.a

lib/libPrg.a: $(GRAPHICS)
	ar r lib/libPrg.a $?
	ranlib lib/libPrg.a

rg/P_graph1.o: rg/P_graph1.c h/P_graphic.h

rg/P_graph2.o: rg/P_graph2.c h/P_graphic.h

rg/P_graph3.o: rg/P_graph3.c h/P_graphic.h

rg/P_graph4.o: rg/P_graph4.c h/P_graphic.h

#
#   And now all together
#

lib/libPILAR.a: $(SPECIFICS) $(GENERICS) P_db.o $(DBFORMATS) \
		$(GRAPHICS)
	ar r lib/libPILAR.a $?
	ranlib lib/libPILAR.a

#
#   Maintenance procedures
#

clean:
	-rm *.o
	-rm db/*.o
	-rm rg/*.o
	-rm core
	-rm mon.out
	-rm *.Addrs
	-rm *.Counts
	-rm *.pixie

clobber: clean
	-rm lib/libPILAR.a
	-rm lib/libPlib.a
	-rm lib/libPgen.a
	-rm lib/libPdb.a
	-rm lib/libPrg.a
	-rm Pilar_?
	-rm Pilar_??
