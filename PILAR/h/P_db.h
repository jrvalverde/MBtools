/**
 *  \file   P_db.h
 *  \brief  Generic database access functions
 *
 *	Header file with the definitions pertaining to module
 *  P_db.c. It is usually #included via P_extern.h
 *
 *	Requires PORTABLE.H for definition of type status.
 *
 *  \see    portable.h
 *
 *  \author     J. R. Valverde	
 *  \date   	8-apr-1990
 */

/** typedef void *db_t;
 *  \brief     Type definition for a generic database 
 *
 *  This typedef hides the actual structure definition of db_t from 
 *  module users. It forces it to be used as a black box.
 */
typedef void *db_t;

/*
 *  DB access functions
 */

/**
 *  \fn     db_t db_open(char *, int);
 *  \brief  Opens a protein database.
 *
 *	It first tries to open the file. If succeeded then it reads
 *   the first line in the database and sets the pointer into
 *   the file appropriately. Also saves the current database name
 *   for possible use later.
 *
 *  \param  db_filnam	The database file name.
 *  \param  db_type 	The type of database (if known) or 0 if unkown.
 *
 *  \return a pointer to a database structure.
 *
 *  \author J. R. Valverde
 */
extern db_t db_open(char *db_filnam, int db_type);

/**
 *  \fn status db_next_protein(db_t db, protein_t *prot);
 *  \brief  Get next protein from the database.
 *
 *	Gets next protein from the database. It sequentially access the
 *  database and reads the protein name and sequence. Also it gets
 *  the index of the protein into the database and obtains the access
 *  keywords.
 *
 *  \param  db	Database to get the protein from.
 *  \param  prot    Structure of type protein to store the read sequence
 *
 *  \return TRUE if we could read a new protein, FALSE otherwise.
 *
 *  \author J. R. Valverde
 */
extern status db_next_protein(db_t, protein_t *);

/**
 *  \fn status db_seek_protein(db_t db, index_t *p_index, protein_t *the_prot);
 *  \brief Get the information about a protein from the database
 *
 *  Gets the information about a protein from the current database.
 *  It access the file directly by using the supplied index.
 *  Obtains all the info about the protein and updates its structure
 *  appropriately
 *
 *  \param  db	Database to get the protein from.
 *  \param  p_index Index of the sequence to retrieve
 *  \param  the_prot	Protein structure to receive the sequence
 *
 *  \return TRUE if success, FALSE if not.
 *
 *  \author J. R. Valverde 
 */
extern status db_seek_protein(db_t, index_t *, protein_t *);

/**
 *  \fn     void db_name(db_t, char *);
 *  \brief  Return the database file name
 *
 *	Returns the current database file name by copying it into
 *  the supplied string.
 *
 *  \param  db	The database whose filename we want to know.
 *  \param  filnam  String to which the name must be copied.
 *
 *  \warning	The routine ASSUMES that this string is long enough.
 *
 *  \author J. R. Valverde
 */
extern void db_name(db_t db, char *filnam);

/**
 *  \fn void db_close(db_t db);
 *  \brief  Close the specified database.
 *
 *  \param  db	Database to be closed.
 *
 *  \author J. R. Valverde
 */
extern void db_close(db_t db);

/**
 *  \fn int db_type(db_t db);
 *  \brief  Return the format of the specified database
 *
 *  \param db	Database to check.
 *
 *  \return The database format
 *
 *  \author J. R. Valverde
 */
extern int db_type(db_t db);

/**
 *  \fn status idx_open(char *file);
 *  \brief Open index file.
 *
 *  \param  file Name of index file.
 *
 *  \return TRUE if success, FALSE if not.
 *
 *  \author J. R. Valverde
 */
extern status idx_open(char *file);

/**
 *  \fn status idx_close();
 *  \brief Close the index file.
 *
 *  \param  file Name of index file.
 *
 *  \return TRUE if success, FALSE if not.
 *
 *  \author J. R. Valverde
 */
extern status idx_close();

/**
 *  \fn status idx_name(char * name);
 *  \brief Return name of the last index file read.
 *
 *  \param  file String to store the name in.
 *
 *  \return TRUE if success, FALSE if not.
 *
 *  \warning	It is assumed that the string passed as argument is long enough.
 *
 *  \author J. R. Valverde
 */
extern status idx_name(char *);

/**
 *  \fn status get_index(index_t *idx);
 *  \brief Read index from index file.
 *
 *  \param  idx Structure to hold the read index.
 *
 *  \return TRUE if success, FALSE if not.
 *
 *  \author J. R. Valverde
 */
extern status get_index(index_t *idx);

/*
 *  Bank is closed now
 */

