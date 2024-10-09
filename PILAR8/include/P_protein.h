/*
 *	P_PROTEIN.H
 *
 *		Header file with the definitions pertaining to module
 *	P_Protein.C. It is usually #included via P_Extern.H
 *
 *		Requires PORTABLE.H for definition of type status.
 *
 *		Designed by:
 *			J. R. Valverde		8  - apr - 1990
 */


extern status read_protein(char *, protein_t *);

extern status get_protein(protein_t *);

extern status protein_fname(char *);

extern char *protein_name(protein_t *);

extern int protein_seq_length(protein_t *);

extern char *protein_sequence(protein_t *);

extern float *protein_feat_values(protein_t *, int);

extern list_t protein_feat_regions(protein_t *, int);

extern void init_protein(protein_t *);

/*
 *	And that's all.
 */
