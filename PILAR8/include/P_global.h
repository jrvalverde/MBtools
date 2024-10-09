/**
 *	\file P_global.h
 *  	\brief Global variables for PILAR	
 *								
 *	SYNOPSIS						
 *		Defines the following globals:			
 *		aa_1[]	Table of one letter code aminoacids	
 *		aa_3[]	Table of three letter code aminoacids	
 *		aa_number	Number of aminoacids.		
 *		aa_table[] Table for finding index into aa_?[]	
 *								
 *	DESCRIPTION						
 *		Module containing global definitions for	
 *	PILAR (Protein Identification and Library Access	
 *	Resource). To include ONLY ONCE in the definitive	
 *	program.						
 *								
 *	CAVEATS							
 *								
 *	\note							
 *		Some compilers will allow to make a module out	
 *	of this file. Others will require to use it as an	
 *	include file to be used ONCE on the program (usually	
 *	at the main module). I choose the last approach as the	
 *	most portable, though the first is more elegant. Check	
 *	your particular system to decide which one to use.	
 *								
 *	SEE ALSO						
 *								
 *	HISTORY							
 *		Creation 27 - jun - 1989 (J. R. Valverde)	
 *
 *  	\author J. R. Valverde								
 */

char aa_1[] = "ACDEFGHIKLMNPQRSTVWXY";

int aa_number = 21;

char *aa_3[] = {"Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His",
		"Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln",
		"Arg", "Ser", "Thr", "Val", "Trp", "***", "Tyr"};

int aa_table[] = {
	-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0,-1, 0, 0, 0,
	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	-1, 1, 3, 2, 3, 4, 5, 6, 7, 8, 0, 9,10,11,12, 0,
	13,14,15,16,17, 0,18,19,20,21, 4, 0, 0, 0, 0, 0,
	 0, 1, 3, 2, 3, 4, 5, 6, 7, 8, 0, 9,10,11,12, 0,
	13,14,15,16,17, 0,18,19,20,21, 4, 0, 0, 0, 0, 0,
};

