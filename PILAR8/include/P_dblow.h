/**
 *  \file P_dblow.h
 *  \brief  Type definitions common to all low level database access operations.
 *
 *  \warning
 *  The included type definitions are made only for use by low 
 *  level database access modules. It should never be used nor
 *  seen on any other modules in the program. ALL DB ACCESS
 *  MUST BE DONE THROUGH THE FUNCTIONS PROVIDED BY P_db.c
 *
 *  \see P_db.c 
 *  \see P_db.h
 *
 *  \author J. R. Valverde	    JRAMON@emduam51.BitNet
 */

/**
 *  \internal
 *  \typedef struct db_st *db_t;
 *  \brief  Type of low-level databases.
 *
 *  \struct db_st
 *  \brief  Internal structure to hold low-level database information.
 */
typedef struct db_st {
    FILE *db_fp;    	    	/**< Pointer to database file */
    int  db_fd;     	    	/**< Database fie descriptor */
    char db_name[MAXPATH];  	/**< Database file name */
    char db_buffer[MAXLINE];	/**< Internal read buffer */
    int  db_linelen;	    	/**< Length of last line read */
    long db_pos;    	    	/**< Current position in file */
    int  db_type;   	    	/**< Type of database */
    status (*db_next_protein)();    /**< Database-specific fecth routine */
} *db_t;

