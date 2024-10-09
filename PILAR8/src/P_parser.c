/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*							    *
*	MODULE:						    *
*	    P_PARSER.C					    *
*							    *
*	SYNOPSIS:					    *
*	    This module contains the parser for the rules   *
*	to be used in the PILAR.			    *
*							    *
*	DESCRIPTION:					    *
*	    The module contains the following entry points: *
*		parse()					    *
*		store_token()				    *
*		retrieve_token()			    *
*							    *
*	CAVEATS:					    *
*							    *
*	NOTES:						    *
*	    Requires the position pointer in the file to be *
*	parsed to be located at the beginning of the rules. *
*	    The tokens read are stored in a linked list	    *
*	that should be accessed via store_token(),	    *
*	retrieve_token() from the outside. Token values are *
*	stored as asciiz strings and SHOULD BE DECONVERTED  *
*	by the caller.					    *
*							    *
*	FILES:						    *
*	    The module depends on the following files:	    *
*		STDIO.H					    *
*		PORTABLE.H				    *
*		P_DEFS.H				    *
*		P_TYPES.H				    *
*		P_EXTERN.H				    *
*		P_LIST.C				    *
*							    *
*	SEE ALSO:					    *
*	    Module P_GETINF.C contains the function that    *
*	will use parse(): read_descriptor().		    *
*							    *
*	DESIGNED BY:					    *
*	    José Ramón Valverde Carrillo. (1989)	    *
*							    *
*	HISTORY:					    *
*	    16 - sep - 1989 J. R. Valverde (implementation  *
*		of basic skeleton)			    *
*	    18 - sep - 1989 J. R. Valverde (It does not use *
*		the information retrieved !)		    *
*	    20 - sep - 1989 J. R. Valverde (error message   *
*		reporting in explanatory way)		    *
*	    22 - sep - 1989 J. R. Valverde (general parse   *
*		function implemented)			    *
*	    24 - sep - 1989 J. R. Valverde (use tokens list *
*		to pass the data. Debugged, seems OK.)	    *
*	     1 - aug - 1991 J. R. Valverde (ported to	    *
*		Ultrix)					    *
*	    19 - dec - 1992 J. R. Valverde (added comment   *
*		support)				    *
*							    *
*	COPYRIGHT:					    *
*	    This is public domain software. If you want to  *
*	use it for commercial purposes you should contact   *
*	the designer. $$ will be appreciated (or Pounds or  *
*	Francs or Marcs or Pts... )			    *
*					    YoEgo.	    *
*							    *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#ifndef DEBUG
#define DEBUG	    0
/* #define DEBUG    6 */
#endif

/*
 *	Tout homme tend à devenir machine. Habitude, méthode, maîtrise,
 *  en fin cela veut dire machine...
 *
 *	En verité c'est la machine qui tend à devenir animal...
 *
 *	Paul Valéry. Cahiers. La Pléïade, 885, 139
 */

/*
 *	GENERAL INCLUDE FILES
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>

/*
 *	PARTICULAR INCLUDE FILES
 */

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"

#include "P_list.h"

/* Uncomment following line in System III, System V and Xenix 3.0 */
/* #include "P_utils.h"		/ * for index */

/*
 *	PRIVATE DEFINITIONS
 */

#define NO_RULES	3
#define MAX_PAT_SIZE	MAXLINE
#define MAX_TOKEN	MAXLINE

static bool is_token(char *token, char *name);
static int get_next_token(int type);
static void get_delimiter(int type);
static void get_value();
static void get_sequence();
static void get_next_word();
static void parser_error(char *msg);


/*
*		Syntax:
*	    ===============
*
*   <rule> ::= { distance | sequence | profile }
*
*   <distance> ::= <dis_token> <pos> - <pos> = <range>
*   <sequence> ::= <seq_token> <seq> <in_token> <pos_range>
*   <profile>  ::= <pro_token> <region> - <values> < <value>
*
*   <dis_token> ::= dis | dist | dista | distan | distanc | distance
*   <seq_token> ::= seq | sequ | seque | sequen | sequenc | sequence
*   <pro_token> ::= pro | prof | profi | profil | profile
*   
*   <pos> ::= ( <value>, <value>, <value> )
*   <range> ::= <value> .. <value>
*   <seq> ::= " <ASCII> "
*   <in_token> ::= in
*   <pos_range> ::= ( <value>, <value>, <range> )
*   <region> ::= ( <value>, <value> )
*   <values> ::= ( <value> [, <value>] )
*   <value>  ::= <float>
*
*   ;<anything> ::= commentary line. Ignored
*   Blank lines are ignored
*   Heading white space is ignored
*   Anything trailing a complete rule is ignored
*   A position of zero means the end of the corresponding region
*
*/

enum {
    DIS_TOKEN,
    SEQ_TOKEN,
    PRO_TOKEN
} RULE_TOKENS;

/*
 *	VARIABLES PRIVATE TO THE MODULE
 */

static char buff[MAX_PAT_SIZE + 1];
static char *pat = buff;
static char token[MAX_TOKEN + 1];
static char *p_token = token;
static FILE *pat_file;
static descriptor_t *d;
static list_t token_list;
static char *rules_table[] = {
    "DISTANCE",
    "SEQUENCE",
    "PROFILE" };


/*
 * FORWARD DEFINITIONS. These are functions private to THIS module 
 */ 

extern void get_integer();
extern void get_delimiter();
extern void get_sequence();
extern void get_value();
extern void parser_error();
extern void get_next_word();

/*----------------------------------------------------------------------
    PARSE
	Parse a file containing the rules to evaluate in PILAR
	The place pointer in the file must be set to the position
	where the rules begin.

    INPUT:  PATTERN_FILE -- File containing the rules
	    DATA -- List into which rules are to be stored
    OUTPUT: NONE

    J. R. VALVERDE
----------------------------------------------------------------------*/

void parse(pattern_file)
FILE *pattern_file;
{
    int result;
    extern FILE *pat_file;

    if (pattern_file == NULL) return;
    pat_file = pattern_file;

    token_list = new_list(); 

    /* Get all information from the file */ 
    do {
	result = get_next_token( NEW_RULE );
    } while (result != EOF);

}

/*--------------------------------------------------------------------
    STORE_TOKEN
	Stores a token in a tokens list inner to the module. The access
    to this list is limited to the use of this function and of
    retrieve_token(). Each item in the list is token. Each token has a
    type and associated data.

    INPUT:  TYPE -- Token type
	    DATA -- String of associated data building up the token.
    OUTPUT: NONE

    J. R. VALVERDE
---------------------------------------------------------------------*/

void store_token(type, data)
flag type;
char *data;
{
    token_t *new_token;
    char *tk_text;

    new_token = (token_t *) malloc( sizeof( token_t ) );
    if (new_token == NULL)
	parser_error("not enough memory");
    if (data == NULL)
	new_token->tk_data = NULL;
    else {
	tk_text = calloc(strlen(data) + 1, sizeof(char));
	if (tk_text == NULL)
	    parser_error("not enough memory");
	strcpy(tk_text, data);
	new_token->tk_data = tk_text;
    }
    new_token->tk_type = type;
    list_append(new_token, token_list);
}

/*--------------------------------------------------------------------
    RETRIEVE_TOKEN
	This function deletes a node from the list of tokens and
    returns the type of the token and its associated data. If the
    list is empty then it frees the list of tokens.

    INPUT:  TYPE -- A pointer to a flag that will get, on output,
	    the type of token retrieved.
    OUTPUT: A pointer to a string containing the token's data
	    TYPE will be NO_MORE_TOKENS if the end of the list has been
	    reached.

    J. R. VALVERDE
--------------------------------------------------------------------*/

char *retrieve_token(type)
flag *type;
{
    token_t *next_token;
    char *tk_text;
    
    if (token_list == NULL) {
	*type = NO_MORE_TOKENS;
	return NULL;
    }
#ifdef V_3
    if (token_list->l_length == 0) {
#else
    if (list_length( token_list ) == 0) {
#endif
	free( token_list );
	token_list = NULL;
	*type = NO_MORE_TOKENS;
	return NULL;
    }
    next_token = (token_t *) list_del_first( token_list );
    *type = next_token->tk_type;
    tk_text = next_token->tk_data;
    free( next_token );
    return ( tk_text );
}

/*---------------------------------------------------------------------

		    STATIC (PRIVATE) FUNCTIONS

---------------------------------------------------------------------*/

static int new_line()
/*
*   gets a new line from the file containing the information about
*   the patterns to be matched.
*	Arguments: none
*	Returns: EOF if end of file, otherwise, the number of characters read
*
*	J. R. Valverde
*/
{
    extern FILE *pat_file;

    if (feof(pat_file))
	return EOF;

    /* ignore blank lines and lines beginning with ';' */
    do {
	if ( fgets(buff, MAX_PAT_SIZE, pat_file) == NULL )
	    return EOF;
	pat = buff;
	while (isspace(*pat)) pat++;
    } while ((*pat == '\0') || (*pat == ';'));

    /* convert string to uppercase */
    for (pat = buff; *pat != '\0'; ++pat)
	if (isalpha(*pat) && islower(*pat))
	    *pat = toupper(*pat);

    /* eliminate \n */
    pat = index(buff, '\n');
    if (pat != NULL) *pat = '\0';

    /* reset pat */
    pat = buff;

    return (strlen(pat));
}

static void skip_blanks()
/*
*   Skips white space trailing from the current position of the pattern
*   Arguments: none
*   Returns:   none
*
*   J. R. Valverde
*/
{
    extern char *pat;

    while (isspace(*pat)) pat++;
}

static int lookup(token, table, no_entries)
char *token;
char *table[];
int no_entries;
/*
*   Looks for the presence of a given token in a specified table.
*   Allows for incomplete token names. In case of ambiguous token
*   specification, it will return the FIRST matching token type.
*   Arguments:	TOKEN -- The putative token
*		TABLE -- The table of admissible tokens
*		NO_ENTRIES -- The size of the table
*   Returns: number of the token in the table, -1 if not found
*
*   J. R. Valverde
*/
{
    int i;
    extern bool is_token();

    if (*token == '\0') return (-1);

    for (i = 0; i < no_entries; ++i)
	if (is_token(token, table[i]))
	    return (i);

    return (-1);
}

static bool is_token(token, name)
char *token, *name;
/*
*   Compares a putative token with a name. Allows for partial
*   specification of the name.
*   Arguments:	TOKEN -- The putative token
*		NAME -- The admissible full name
*   Returns: TRUE if the token matches the name, FALSE if not
*
*   J. R. Valverde
*/
{
    while (*token++ == *name++);
    if (*token == '\0')
	return TRUE;
    else
	return FALSE;
}

static void init_token()
{
    int i;

    p_token = token;
    for (i = 0; i <= MAX_TOKEN; ++i)
	token[i] = '\0';
}

static int get_next_token(type)
int type;
/*
*   Tries to get next token, which must be of the specified type.
*   If it can not be found then the program aborts after printing
*   an error message.
*   Arguments:	TYPE -- the type of token to get
*   Returns: TYPE if all has gone OK.
*
*   J. R. Valverde
*/
{
    extern void init_token();
    extern int new_line();
    extern void skip_blanks();

    init_token();

    switch (type) {
	case NEW_RULE:
#if DEBUG == 6
	    fprintf(stderr, "NEW_RULE:\n");
#endif
	    if (new_line() == EOF) return ( EOF );
	    get_next_word();
	    store_token(type, NULL);
	    switch (lookup(token, rules_table, NO_RULES)) {
		case (DIS_TOKEN):
		    get_next_token( DISTANCE ); break;
		case (SEQ_TOKEN):
		    get_next_token( SEQUENCE ); break;
		case (PRO_TOKEN):
		    get_next_token( PROFILE );  break;
      		default:
		    parser_error("Unknown rule"); break;
	    }
 	    break;

	case DISTANCE:
#if DEBUG == 6
	    fprintf(stderr, "	DISTANCE\n");
#endif
	    store_token(type, NULL);
	    get_next_token( POSITION );
	    get_next_token( HYPHEN );
	    get_next_token( POSITION );
	    get_next_token( EQUAL );
	    get_next_token( RANGE );
	    break;

      	case SEQUENCE:
#if DEBUG == 6
      	    fprintf(stderr, "	SEQUENCE\n");
#endif
	    store_token(type, NULL);
	    get_next_token( SEQ );
	    get_next_token( IN_TOKEN );
	    get_next_token( POS_RANGE );
      	    break;

       	case PROFILE:
#if DEBUG == 6
	    fprintf(stderr, "	PROFILE\n");
#endif
	    store_token(type, NULL);
	    get_next_token( REGION );
	    get_next_token( HYPHEN );
	    get_next_token( VALUES );
	    get_next_token( LESSTHAN );
	    get_next_token( VALUE );
	    break;

	case POSITION:
#if DEBUG == 6
	    fprintf(stderr, "\t\tposition:\n");
#endif
	    get_next_token( LPAREN );
	    get_next_token( VALUE );
	    get_next_token( COLON );
	    get_next_token( VALUE );
	    get_next_token( COLON );
	    get_next_token( VALUE );
	    get_next_token( RPAREN );
	    break;

	case RANGE:
#if DEBUG == 6
	    fprintf(stderr, "\t\trange\n");
#endif
	    get_next_token( VALUE );
	    get_next_token( POINTS );
	    get_next_token( VALUE );
	    break;

	case SEQ:
#if DEBUG == 6
	    fprintf(stderr, "\t\tsequence\n");
#endif
	    get_sequence();
	    store_token( type, token );
	    break;

	case IN_TOKEN:
#if DEBUG == 6
	    fprintf(stderr, "\t\tin\n");
#endif
	    get_next_word();
	    if (strcmp(token, "IN") != 0)
		parser_error("IN expected");
	    break;

      	case POS_RANGE:
#if DEBUG == 6
	    fprintf(stderr, "\t\tpos_range\n");
#endif
	    get_next_token( LPAREN );
	    get_next_token( VALUE );
	    get_next_token( COLON );
	    get_next_token( VALUE );
	    get_next_token( COLON );
	    get_next_token( RANGE );
	    get_next_token( RPAREN );
	    break;

      	case REGION:
#if DEBUG == 6
	    fprintf(stderr, "\t\tregion\n");
#endif
	    get_next_token( LPAREN );
	    get_next_token( VALUE );
	    get_next_token( COLON );
	    get_next_token( VALUE );
	    get_next_token( RPAREN );
	    break;

      	case VALUES:
#if DEBUG == 6
	    fprintf(stderr, "\t\tvaules\n");
#endif
	    get_next_token( LPAREN );
	    get_next_token( VALUE );
	    skip_blanks();
	    while (*pat == ',') {
		pat++;
		get_next_token( VALUE );
		skip_blanks();
	    }
	    get_next_token( RPAREN );
	    break;

      	case VALUE:
#if DEBUG == 6
	    fprintf(stderr, "\t\tvalue");
#endif
	    get_value();
	    store_token( type, token );
	    break;

	case HYPHEN:
	case EQUAL:
	case LPAREN:
	case RPAREN:
	case COLON:
	case POINTS:
      	case LESSTHAN:
#if DEBUG == 6
	    fprintf(stderr, "\t\tdelimiter\n");
#endif
	    get_delimiter(type);
	    break;
	default:
	    parser_error("unknown token");
    }

    return (type);
}

static void get_delimiter(type)
int type;
{
    extern char *p_token, *pat;
    static char delimiter[] = "-=(),.<";

    skip_blanks();
    if (*pat != delimiter[ type - HYPHEN ])
	parser_error("bad delimiter");
    token[0] = *pat;
    token[1] = '\0';
    if (type == POINTS) {
	/* we need another point */
	pat++;
	if (*pat != delimiter[ type - HYPHEN ])
	    parser_error("bad delimiter");
	token[1] = *pat;
	token[2] = '\0';
    }
    p_token = token;
    pat++;
}


static void get_value()
/*
*   gets a float value from the input stream
*   Arguments: none
*   Returns: none
*
*   J. R. Valverde
*/
{
    skip_blanks();
    p_token = token;
    /* first check for minus sign */
    if (*pat == '-')
	*p_token++ = *pat++;
    if (! isdigit( *pat ))
	parser_error("numerical value expected");

    do {
	*p_token++ = *pat++;
    } while (isdigit(*pat));

    /* check for fractional part: we use a point to distinguish
    from colons which are recognized as separators. */
    if (*pat == '.') {
     	do {
	    *p_token++ = *pat++;
	} while (isdigit(*pat));
    }

    /* check for exponent */
    if (*pat == 'E') {
	*p_token++ = *pat++;
	/* check for minus sign */
	if (*pat == '-')
	    *p_token++ = *pat++;
	if (! isdigit(*pat))
	    parser_error("missing exponent");
	do {
	    *p_token++ = *pat++;
	} while (isdigit(*pat));
    }
    *p_token = '\0';
    p_token = token;

#if DEBUG == 6
    fprintf(stderr, "\t\t\treal value: %s\n", token);
#endif
}

static void get_sequence()
/*
*   gets an ascii sequence from the input stream
*   Arguments: none
*   Returns: none
*
*   J. R. Valverde
*/
{
    char *tmp;

    skip_blanks();
    if (*pat++ != '"')
	parser_error("quotes expected");
    p_token = token;
    do {
	*p_token++ = *pat++;
    } while (*pat != '\"' && *pat != '\0');
    *p_token = '\0';
    p_token = token;
    pat++;
#if DEBUG == 6
    fprintf(stderr, "\t\t\tsequence: %s\n", token);
#endif
}

static void get_next_word()
/*
*   gets next word from the input stream
*   Arguments: none
*   Returns: none
*
*   J. R. Valverde
*/
{
    skip_blanks();
    if (*pat == '\0')
	parser_error("unexpected end of line");
    p_token = token;
    do {
	*p_token++ = *pat++;
    } while (isalnum(*pat));
    *p_token = '\0';
    p_token = token;
}

static void parser_error(msg)
char *msg;
/*
*   Prints an error message showing the last token read, the
*   actual position of the pointer in the pattern and a possible
*   explanation of the problem. Then it aborts the program.
*   Arguments:	MSG -- explanatory message
*   Returns:	Exit to the operating system with a letal error indication
*
*   J. R. Valverde
*/
{
    int i;

    extern char buff[];
    extern char *pat;
    extern char token[];

    fprintf(stderr, "%s\n", buff);
    for (i = 1; i < (pat - buff); ++i)
	fprintf(stderr, "-");
    fprintf(stderr, "^\n");
    fprintf(stderr, "ERROR: %s -- ' %s '\n", msg, token);
#if DEBUG == 6
    return;
#endif
    exit(LETAL_ERROR);
}

/*--------------------------------------------------------------------
			D E B U G G I N G
--------------------------------------------------------------------*/

#if DEBUG == 6


main()
{
    char name[100];
    FILE *pattern_file;
    flag type;
    char *data;
    int i;

    printf("File to parse: ");
    scanf("%s", name);                     
    pattern_file = fopen(name, "rw+");
/*
    do {
	get_next_token( NEW_RULE );
    }
    while ( ! feof(pattern_file) );
*/
    parse(pattern_file);
    fclose(pattern_file);
#ifdef V_3
    fprintf(stderr, "\nNumber of tokens: %d\n", token_list->l_length);
#else
    fprintf(stderr, "\nNumber of tokens: %d\n", list_length(token_list));
#endif

    i = 0;
    do {
	data = retrieve_token( &type );
	if (data == NULL) fprintf(stderr, "NULL\n");
	else fprintf(stderr, "%s\n", data);
	free(data);
	i++;
    } while (type != NO_MORE_TOKENS);
    fprintf(stderr, "\nNumber of tokens: %d\n", i);
}

#endif
