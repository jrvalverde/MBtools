
/* Part 2 of torture-test of the memory controller */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>

# include "MCon.h"

# define MAXTYPES	20
# define MAXINSTANCES	500
# define BASESIZE	25
# define MAXEXTRA	25
# define NPASSES	500

void torture(void);
void error(char *);
void makeNewType(void);
void makeNewObject(void);
void freeSomeObject(void);
int randle(int);
void debugCon(char *);
struct {
    MCon mCon;
    long theID;
} regType[MAXTYPES];
struct {
    void *mObj;
    long typeID;
} mBag[MAXINSTANCES];
int nInstances = 0;
int nTypes = 0;
long serialID = 1000;
void torture()
{
    counter orphans = 0;
    size_t sizeItem;
    int typeIdx, objIdx;
    int nPotential, nActual;
    long tID;
    int i, x;
    long looper;
    char msg[64];
    for (looper = 0; looper < NPASSES; looper++) {
	if (orphans) error("orphans have been created");
	if (!nTypes) {
	    if (nInstances) error("instances but no active types");
	    makeNewType();
	}
	x = randle(10);
	switch (x) {
	case 0 :		/* make a report */
	    sprintf(msg, "loop %ld : %d instances of %d types\n",
	      looper,
	      nInstances,
	      nTypes);
	    debugCon(msg);
	    break;
	case 1 :	    /* register a new type if room */
	    makeNewType();
	    break;
	case 2 :	    /* make a new object if room */
	    makeNewObject();
	    break;
	case 7 :	    /* make many new objects */
	    nPotential = ((MAXINSTANCES - nInstances) >> 2) + 7;
	    nActual = randle(nPotential);
	    while (nPotential--)
		makeNewObject();
	    break;
	case 3 :	    /* free an object if any exist */
	    freeSomeObject();
	    break;
	case 8 :	    /* free many objects */
	    nPotential = ((MAXINSTANCES - nInstances) >> 3) + 3;
	    nActual = randle(nPotential);
	    while (nPotential--)
		freeSomeObject();
	    break;
	case 4 :      /* purge free list of a type if any exist */
	    if (nTypes) {
		typeIdx = randle(nTypes);
		purgeMCon(regType[typeIdx].mCon);
	    }
	    break;
	case 9 :                /* purge a number of free lists */
	    if (nTypes) {
		nActual = randle(nTypes);
		for (i = 0; i < nActual; i++) {
		    typeIdx = randle(nTypes);
		    purgeMCon(regType[typeIdx].mCon);
		}
	    }
	    break;
	case 5 :     /* free all of a registered type if any exist */
	    if (nTypes) {
		typeIdx = randle(nTypes);
		tID = regType[typeIdx].theID;
		for (i = 0; i < nInstances; )
		    if (mBag[i].typeID == tID) {
			disposeInstance(mBag[i].mObj);
			nInstances--;
			mBag[i] = mBag[nInstances];
		    }
		    else i++;
		    /* and maybe kill the controllor */
		if (randle(13) > 7) {
		 orphans += disposeMCon(regType[typeIdx].mCon);
		 nTypes--;
		 regType[typeIdx] = regType[nTypes];
		}
	    }
	    break;
	case 6 :             /* free all instances of all types */
	    if (randle(20) < 15) break;
	    for (i = 0; i < nInstances; i++)
		disposeInstance(mBag[i].mObj);
	    nInstances = 0;
		    /* kill some controllors */
	    for (i = 0; i < nTypes;)
		if (randle(15) > 13) {
		      orphans += disposeMCon(regType[i].mCon);
		      nTypes--;
		      regType[i] = regType[nTypes];
		}
		else i++;   
	       /* and maybe kill all the rest and shut down! */
	    if (randle(20) > 17) {
		orphans += disposeAllMCons();
		nTypes = 0;
	    }
	    break;
	}
    }
    printf("\n%ld passes...\n", NPASSES);
}
void makeNewType()
{
    MCon tCon;
    size_t sizeItem;
    if (nTypes < MAXTYPES) {
	sizeItem = BASESIZE + randle(MAXEXTRA);
	tCon = newMCon(sizeItem);
	if (!tCon) {
	    error("could not make controllor");
	    return;
	}
	regType[nTypes].mCon = tCon;
	regType[nTypes].theID = serialID++;
	nTypes++;
    }
}
void makeNewObject()
{
    int typeIdx;
    void *tObj;
    if ((nInstances < MAXINSTANCES) && nTypes) {
	typeIdx = randle(nTypes);
	tObj = newInstanceOf(regType[typeIdx].mCon);
	if (!tObj) {
	    error("could not get object memory");
	    return;
	}
	mBag[nInstances].mObj = tObj;
	mBag[nInstances].typeID = regType[typeIdx].theID;
	nInstances++;
    }
}
void freeSomeObject()
{
    int objIdx;
    if (nInstances) {
	objIdx = randle(nInstances);
	disposeInstance(mBag[objIdx].mObj);
	nInstances--;
	mBag[objIdx] = mBag[nInstances];
    }
}
void error(char *s)                    /* print error message and hang */
{
    char hang[12];
    printf("\nERROR : %s\n", s);
    scanf("%s", hang);
}
randle(int n)	                     /* random number up to not including n */
{
    if (n <= 0) error("bogus randle call");
    return (rand() % n);
}


