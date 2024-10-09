/* Header for memory control package 
 C 1989 Robert A. Moeser, all rights reserved */

/* some things probably in your standard headers... */
#ifndef __size_t
#define __size_t
typedef long size_t;
#endif

void *malloc(size_t);
void free(void *);
int printf(char *, ...);
int scanf(char *, ...);
int rand(void);

/* end of things probably in your standard headers... */
typedef long counter;	    /* a bit excessive, perhaps... */
typedef struct _mc {
    struct _mc *next;	/* used for doubly-linked list */
    struct _mc *prev;	/* of all mcs */
    size_t itemSize;	/* the size of the object this MCon controls */
    void *freeList; 	/* a free list of objects */
    counter nGiven; 	/* number of items handed out */
    counter nOut;	/* number still out there somewhere */
    counter nFree;	/* length of free list */
} *MCon;

MCon newMCon(size_t theItemSize);   /* 0 = could not get storage */
void *newInstanceOf(MCon theCon);   /* 0 = could not get storage */
void disposeInstance(void *theItem);
void purgeMCon(MCon theCon);
counter disposeMCon(MCon theCon);    /* number of orphans created by action */
void purgeAllMCons(void);
counter disposeAllMCons(void);	     /* number of orphans created by action */



