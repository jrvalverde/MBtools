/*  A Memory Controller. C 1989 Robert A. Moeser, all rights reserved */

# include "MCon.h"

static void linkIn(MCon theCon);
static void linkOut(MCon theCon);

void debugCon(char *);
void printCon(MCon theCon);

static MCon controlCon = 0;	/* 0 -> not yet initialized */

/* Initialize the memory control "package." Returns 1 for successful 
initialization or already initialized; returns 0 if it fails. Called 
automatically by newMCon, but you can call it explicity if you like. */

int initMControl()
{
    if (controlCon) return (1);	    /* already initialized! */
    controlCon = (MCon) malloc(sizeof *controlCon);
    if (!controlCon) return (0);    /* sad but true! */
    controlCon->itemSize = sizeof *controlCon;
    controlCon->nGiven = controlCon->nOut = controlCon->nFree = 0;
    controlCon->freeList = (void *) 0;
    controlCon->next = controlCon->prev = controlCon;
    return (1);			    /* OK! */
}

/*  Register a new type (object) for management by the memory controller
takes size of the object as an argument and returns an MCon object, which 
is zero in case of failure. The Mcon object is used later to manage creation 
of all instances of the new object, to maintain a free list for the object 
and to keep track of demand for the object */

MCon newMCon(size_t theItemSize)
{
    MCon t;
    int k;
    k = initMControl();
    if (!k) return ((void *) 0);
    t = newInstanceOf(controlCon);
    if (!t) return ((void *) 0);
    t->itemSize = theItemSize;
    t->nGiven = t->nOut = t->nFree = 0;
    t->freeList = (void *) 0;
    linkIn(t);
    return (t);
}

/* Create a new object. Takes the object's MCon as an argument (previously 
created by newMCon). Returns a pointer to a new instance of the object */
void *newInstanceOf(MCon theCon)
{
    void *t;
    if (theCon->freeList) { 	/* any on free list? */
	t = theCon->freeList;
	theCon->freeList = *((MCon *) t);
	theCon->nFree--;
    }
    else	                /* nope, go malloc one */
	t = malloc(theCon->itemSize + (sizeof theCon));
    if (!t) return ((void *) 0);    /* allocator failure... */
    *((MCon *) t) = theCon; 	/* remember the controller */
    theCon->nGiven++;
    theCon->nOut++;
    return ((void *) ((MCon *) t + 1));
}

/* Dispose of an object. Takes a pointer to the object to dispose of
the storage is kept by the object's MCon on a free list for reuse. */
void disposeInstance(void *theItem)
{
    MCon t;
    void *x;
    int g;
    t = *((MCon *) theItem - 1);    /* recover controllor */
    t->nOut--;
    t->nFree++;
    *((MCon *) theItem - 1) = t->freeList;
    t->freeList = (MCon *) theItem - 1;
}

/* Purge the free list for an MCon. Takes the MCon whose free list should be 
returned to the system storage allocator */
void purgeMCon(MCon theCon)
{
    void *bop, *bop2;
    bop = theCon->freeList;
    while (bop) {
	bop2 = *((void **) bop);
	free(bop);
	bop = bop2;
    }
    theCon->freeList = (void *) 0;
    theCon->nFree = 0;
}

/* Unregister a type. Takes the MCon object that is no longer needed
returns the number of "orphans" created. Since active instances are entirely 
the responsibility of clients of the memory control package, disposing of 
an MCon can mean that there may be outstanding objects of the 
no-longer-existing type. This is not strictly an error, but since there is 
now no way to free the storage used by the orphans should be a warning sign. */

counter disposeMCon(MCon theCon)
{
    counter orphans = 0;    /* did this dispose create "orphans" ? */
    orphans = theCon->nOut;
    purgeMCon(theCon);	/* goodbye the free list */
    linkOut(theCon);
    disposeInstance((void *) theCon);
    return (orphans);	/* number of "orphans", should be 0! */
}

/* Purge the free list of every object type known to the memory controller
a client might call this in a desperate attempt to satisfy a memory request */
void purgeAllMCons(void)
{
    MCon bop;
    bop = controlCon;
    do {
	purgeMCon(bop);
	bop = bop->next;
    } while (bop != controlCon);
}

/* Unregister all types. This action can create orphans in the same sense as 
disposeMCon above, and so returns the number. It really ought to be zero 
unless the client has created objects meant to endure until program 
termination. Since all types are deactivated the additional step of 
deactivating the memory controller itself is taken. All storage used by the 
package is returned to the system storage allocator. */

counter disposeAllMCons(void)
{
    counter orphans = 0;
    MCon bop, bop2;
    bop = controlCon->next;
    while (bop != controlCon) {
	bop2 = bop->next;
	orphans += disposeMCon(bop);
	bop = bop2;
    }
    purgeMCon(controlCon);
    free((void *) controlCon);
    controlCon = (void *) 0;	/* de-initialize mc */ 
    return (orphans);	/* total number of "orphans", should be 0! */
}

/* Internal Routines */

/* Link a new MCon into the doubly-linked list of all known MCons
the list head and tail is the MCon for MCons, so there is always at least one
item on the list and the first item is always the MCon's MCon */
static void linkIn(MCon theConToAdd)
{
    theConToAdd->next = controlCon->next;
    theConToAdd->prev = controlCon;
    controlCon->next = theConToAdd;
    (theConToAdd->next)->prev = theConToAdd;
}

/* Unlink an MCon (called upon destruction of an Mcon */
static void linkOut(MCon theConToDel)
{
    (theConToDel->prev)->next = theConToDel->next;
    (theConToDel->next)->prev = theConToDel->prev;
}

/* Debugging and Utility Routines */

/* Print a list of known MCons with statistics. Takes a tag to accompany 
printout as an argument. */
void debugCon(char *s)
{
    MCon bop;
    printf("%s", s);
    if (!controlCon) {
	printf(" -mc is OFF!\n");
	return;
    }
    bop = controlCon;
    do {
	printCon(bop);
	bop = bop->next;
    } while (bop != controlCon);
}

/* Print one MCon. Make a quick check on the integrity of the internal 
data structure. */
void printCon(MCon theCon)
{
    char *freeBop;
    int freeCount;
    int OK;
    freeBop = theCon->freeList;
    freeCount = 0;
    while (freeBop) {
	freeCount++;
	freeBop = *((char **) freeBop);
    }
    OK = freeCount == theCon->nFree;
    printf("%lx = %lu\t%ld\t%ld\t%ld\t%lx\t%s\n",
	theCon,
	theCon->itemSize,	
	theCon->nGiven,
	theCon->nOut,
	theCon->nFree,
	theCon->freeList,
	OK ? "<ok>" : "<NOK>");
}


