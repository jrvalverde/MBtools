/**
 *	\file P_desc.h
 *  	\brief	Interface module for protein descriptors.
 *
 *		Header file with the definitions pertaining to module
 *	P_descr.c. It is usually #included via P_extern.h
 *
 *	\pre	Requires PORTABLE.H for definition of type status.
 *
 *	\see portable.h
 *
 *	\author	J. R. Valverde		
 *  	\date 8-apr-1990
 */

/**
 *  \fn int read_feature(char *file, feature_t *f)
 *  \brief Reads data values from the especified feature values file.
 *
 *  \param file Name of file containing values.
 *  \param  f	Feature table filled with the values.
 *
 *  \note Values must be stored in the following format:
 *	"AminoAcid (1 letter code) HTAB Value CR/LF"
 *
 *  \return Number of amino acid values read.
 */
extern int read_feature(char *file, feature_t *f);

/**
 *  \fn status get_feature(feature_t *f)
 *  \brief Prompt user for feature file and read it.
 *
 *   Asks the user for the feature to be used in the search
 *   and reads it from the file.
 *
 *  \param  f	Table to read data values into.
 *  
 *  \return SUCCESS if read's been OK. FAIL otherwise.
 */
extern status get_feature(feature_t *f);

/**
 *  \fn status feature_name(char *file);
 *  \brief Returns the name of the last feature read.
 *
 *  \param  file String to store the name.
 *
 *  \return SUCCESS if so, FAIL if there is not current feature name.
 *
 *  \warning This function is HIGHLY SYSTEM DEPENDENT.
 *          Supposes that f_name is long enough to hold
 *	the name.
 */
extern status feature_name(char *file);

/**
 *  \fn status read_descriptor(char *file, descriptor_t *d)
 *  \brief  Read descriptor
 *
 *  \param file Name of file containing the descriptor definition.
 *  \param d	Structure to store the descriptor in.
 *
 *  \return TRUE if all went well, FALSE otherwise.
 */
extern status read_descriptor(char *file, descriptor_t *d);

/*
 *	And that's all.
 */
