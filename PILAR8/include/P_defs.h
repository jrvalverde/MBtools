/**								
 *   \file  P_defs.h						
 *   \brief Constants used in PILAR
 *						
 *   Constant definitions for the Protein Identification and	
 *   Library Access Resource.					
 *								
 *								
 *  \author José Ramón Valverde Carrillo.				
 *								
 *   HISTORY:							
 *	J. R. Valverde (1989):					
 *	    Updated as needed while developing the project.	
 *	J. R. Valverde ( 15 - aug - 1989):			
 *	    Version 1 concluded.				
 *	J. R. Valverde ( 20 - aug - 1989):			
 *	    Version 2 terminated.				
 *	J. R. Valverde ( 08 - oct - 1989):			
 *	    Version 3 finished.	
 *
 *  \nosubgrouping				
 */

/**
 *  \name General
 *   General interest definitions.  
 */

/**@{*/

/** Longest input line			*/
#define	MAXLINE	    256		

/** Number of index keywords		*/
#define	KEYS_NUMBER 1		

/** Length of keyword			*/
#define KEYLEN	    10		

/** Maximum protein name length		*/
#define	PNAME_LEN   80		

/** Maximum sequence length		*/
#define	MAX_SEQ_LEN 10000	

/** Maximum number of regions 		*/
#define MAX_REGIONS 100		

/** Maximum number of value categories	*/
#define MAX_CLASSES 10		

/** Number of aminoacids			*/
#define AA_NUMBER   21		

/** The region has already been used	*/
#define ASSIGNED    1		

/** This rule has already been checked	*/
#define EXPLORED    1		


#ifndef OK
/** The protein passes this rule		*/
#define	OK	    1		
#endif

/** Cannot align a region with a pattern */
#define NOT_ALIGNED -1		

/** maximum float			*/
#define F_INFINITY  9.99e+36	

/** A macro to safely allocate space for \a x	*/
#define NEW(x)	    getmem( sizeof( x ) )

/** file extension for protein files	*/
#define P_EXT	    ".aa"	

/** file extension for descriptor files	*/
#define D_EXT	    ".dxr"	

/** default file extension for feat file    */
#define F_EXT	    ".nor"	

/**@}*/

/**
 *  \name db
 * Definitions of low level database access types.
 *
 *  \see P_db.h
 *  \see P_db.c
 */

/**@{ */
#define DB_UNKNOWN	    0
#define DB_STANFORD	    1
#define DB_INTELLIGENETICS  2
#define DB_STRIDER	    3
#define DB_GENBANK	    4
#define DB_NBRF		    5
#define DB_PIR		    6
#define DB_PEARSON	    7
#define DB_EMBL		    8
#define DB_SWISSPROT	    9
#define DB_UWGCG	    10
#define DB_ZUKER	    11
#define DB_FITCH	    12
#define DB_PREER	    13

/**@}*/

/**
 *  \name Queue
 * Definitions for queue manipulation.
 *
 *  \see P_queue.h
 *  \see P_queue.c
 */

/**@{*/ 
                   
/** size of static circular queue */
#define QUEUE_SIZE  128		

/** error: queue is full	    */
#define Q_FULL	    -1		

/** error: queue is empty    */
#define Q_EMPTY	    NULL	

/**@}*/

/**
 * \name Token
 * Definitions for token manipulation while parsing files .
 *
 * \see P_parser.h
 * \see P_parser.c
 */

/**@{*/

/**
 * \enum TOKENS
 * \brief The kinds of tokens known to the parser.
 */
enum {
    NEW_RULE,
    DISTANCE,
    SEQUENCE,
    PROFILE,
    POSITION,
    RANGE,
    SEQ,
    IN_TOKEN,
    POS_RANGE,
    REGION,
    VALUES,
    VALUE,
    HYPHEN,
    EQUAL,
    LPAREN,
    RPAREN,
    COLON,
    POINTS,
    LESSTHAN
}   TOKENS;

/** No more tokens remain to be processed */
#define NO_MORE_TOKENS	-1  

/**@}*/

/**
 * \name Exit
 * Definitions for program termination 
 */

/**@{*/

#define BYE 0

#ifdef VMS
/** 
 * \note The following definitions should be modified according to the
 *  system being used:
 *
 *  	\li EXIT_OK should be 0 for UNIX, 1 for VAX-VMS.
 *  	\li MAXPATH should be set to the appropriate value for the corresponding
 *	system.
 */

/**
 *  \note Vax-VMS Condition codes (exception codes, not only exit codes):
 *
 *  \code
 *	31   28 27                16 15                 3 2   0
 *	+------+--------------------+--------------------+----+
 *	|      |                    |                    |    |
 *	+------+--------------------+--------------------+----+
 *	Control Facility number	     Message number	  Severity
 *  \endcode
 *
 *  	Severity - Severity of the error condition. Bit 0 indicates
 *	success when set. Bits 1 & 2 distinguish degrees of success
 *	or failure:
 *	    \li 0	Warning.    Execution continues. Unpredictable results
 *	    \li 1	Success.    Execution continues. Expected results
 *	    \li 2	Error.	    Execution continues. Erroneous results
 *	    \li 3	Informat.   Execution continues. Information message
 *	    \li 4	Fatal.	    Execution terminates. No output
 *	    \li 5 - 7	    Reserved for DEC
 *
 *  	Message number - Number of the message associated with the
 *	error condition (the msg. may or not be displayed).
 *
 *	Facility number - Number identifying the program originating
 *	the message. Bit 27 is set for user programs and clear for DEC's
 *
 *	Control - Bit 28 inhibits the display of the error message. Bits
 *	29 - 31 are reserved for DEC.
 *
 * \li BIT_0    0x01
 * \li BIT_1    0x02
 * \li BIT_2    0x04
 * \li BIT_15   0x8000
 * \li BIT_16   0x10000
 * \li BIT_27   0x8000000
 * \li BIT_28   0x10000000
 *
 *		SevCode/User progr./Don't display
 */

#define VMS_WARNING (0 | 0x8000000 | 0x10000000)
#define VMS_SUCCESS (1 | 0x8000000 | 0x10000000)
#define VMS_ERROR   (2 | 0x8000000 | 0x10000000)
#define VMS_INFO    (3 | 0x8000000 | 0x10000000)
#define VMS_FATAL   (4 | 0x8000000 | 0x10000000)

/** Successful termination	*/
#define EXIT_OK	    VMS_SUCCESS	    

/** Informative message	*/
#define BAD_ARGS    VMS_INFO	    

/** Letal error message	*/
#define LETAL_ERROR VMS_FATAL	    

/** Max length of file name	*/
#define	MAXPATH	    128

/** Options delimiter char	*/
#define OPTION_DELIMITER    '/'	    

#define DRIVECHAR	    ':'

/** use \ for MSDOS, : for Mac */
#define BEGDIRCHAR	    '['	    

#define DIRDELIMITER	    '.'

#define ENDDIRCHAR	    ']'

#define FILEXTCHAR	    '.'


#else
#ifdef UNIX

#define UNIX_SUCCESS	0

#define EXIT_OK	    UNIX_SUCCESS
#define BAD_ARGS    1
#define LETAL_ERROR 2
#define MAXPATH	    128
#define OPTION_DELIMITER    '-'
#define DRIVECHAR	    ':'
#define BEGDIRCHAR	    '/'
#define DIRDELIMITER	    '/'
#define ENDDIRCHAR	    '/'
#define FILEXTCHAR	    '.'

#endif
#endif
                               
/**@}*/
