/**
 *	\file P_extern.h
 *  	\brief external definitions for PILAR	
 *								
 *	SYNOPSIS						
 *								
 *	DESCRIPTION						
 *	    Contains all the external definitions for the	
 *	different modules of the Protein Identification and	
 *	Library Access Resource. Include it for making all	
 *	the needed functions globally available.		
 *								
 *	CAVEATS							
 *	    This is dynamically modified as the project is	
 *	being built up.						
 *								
 *	FILES							
 *	    PORTABLE.H						
 *	    P_GLOBALS.C						
 *	    P_PROTEIN.C						
 *	    P_DESCRIPTOR.C					
 *	    P_PARSER.C						
 *	    P_DBANK.C						
 *	    P_INIT.C						
 *	    P_REGION.C						
 *	    P_VALUES.C						
 *	    P_UTILS.C						
 *	    P_LIST.C						
 *	    P_DLIST.C						
 *	    P_STACK.C						
 *	    P_QUEUE.C						
 *	    P_SORT.C						
 *			       					
 *	\note							
 *	    Some definitions may be never be used. They are 	
 *	there just for historical reasons.			
 *								
 *	\see						
 *	    It may be worth having a look to the mentioned	
 *	files to get more information about these routines.	
 *								
 *	\author						
 *	    José Ramón Valverde Carrillo.			
 *								
 *	HISTORY							
 *	    8  - jul - 1989 J. R. Valverde			
 *	    26 - sep - 1989 J. R. Valverde			
 *	    8  - oct - 1989 J. R. Valverde			
 *	    8  - apr - 1990 J. R. Valverde			
 *								
 *	COPYRIGHT						
 *	    As far as I know this is supposed to be public	
 *	domain software. Please, do not use it for commercial	
 *	purposes without the written permission of the author.	
 *	(I could become very angry if I became very hungry just	
 *	because of you).					
 *						    J. R.	
 *
 */

/** P_global.c
 *  Global variables from P_GLOBAL.C
 *
 *	Some compilers do not allow a module made only of a
 *  data segment without code segment. It is for this reason
 *  that module P_GLOBAL.C is actually P_GLOBAL.H and must
 *  included once at the main module. Therefore, I can't build
 *  a header file P_GLOBAL.H, and these definitions come in
 *  this file, P_EXTERN.H
 */

/** */
extern char aa_1[];

extern int aa_number;

extern char *aa_3[];

extern int aa_table[];

/*
 *  Routines from P_PROTEIN.C
 */

#ifndef _P_PROTEIN_
#include "P_protein.h"
#endif

/*
 *  Routines from P_DESCRIPTOR.C
 */

#ifndef _P_DESCRIPTOR_
#include "P_desc.h"
#endif

/*
 *  Routines from P_PARSER.C
 */

#ifndef _P_PARSER_
#include "P_parser.h"
#endif

/*
 *  Routines from P_DBANK.C
 */

/* D.Bank: Visa, American Express and MasterCard accepted */

#ifndef _P_DBANK_
#include "P_db.h"
#endif

/*
 *  Routines from P_INIT.C
 */

#ifndef _P_INIT_
#include "P_init.h"
#endif

/*                 
 *  Routines from P_REGION.C
 */

#ifndef _P_REGION_
#include "P_region.h"
#endif

/*
 *  Routines from P_VALUES.C
 */

#ifndef _P_VALUES_
#include "P_values.h"
#endif

/*
 *  Routines from P_THINK.C
 */

/*  Don't worry be happy ! */

#ifndef _P_THINK_
#include "P_think.h"
#endif

/*
 *  Routines from P_UTILS.C
 */

#ifndef _P_UTILS_
#include "P_utils.h"
#endif

/*
 *  Routines from P_LIST.C
 */

#ifndef _P_LIST_
#include "P_list.h"
#endif

/*
 *  Routines from P_DLIST.C
 */

#ifndef _P_DLIST_
#include "P_dlist.h"
#endif

/*
 *  Routines from P_STACK.C
 */

#ifndef _P_STACK_
#include "P_stack.h"
#endif

/*
 *  Routines from P_QUEUE.C
 */

#ifndef _P_QUEUE_
#include "P_queue.h"
#endif

/*
 *  Routines from P_GTREE.C
 */

#ifndef _P_GTREE_
#include "P_gtree.h"
#endif

/*
 *  Routines from P_STREE.C
 */

#ifndef _P_STREE_
#include "P_stree.h"
#endif

/*
 *  Routines from P_AVL.C
 */

#ifndef _P_AVL_
#include "P_AVL.h"
#endif

/*
 *  Routines from P_SORT.C
 */

#ifndef _P_SORT_
#include "P_sort.h"
#endif


/*
 *	END OF MODULE P_EXTERN.H
 */

