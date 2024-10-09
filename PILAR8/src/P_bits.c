/**
 *	\file P_bits.c
 *  	\brief	Bit manipulation routines.
 *
 *	Módulo con rutinas de manipulación de vectores de bits
 *	extraído de una librería de rutinas para CP/M
 *
 *	CONTENIDO:
 *		bv_new(bv, bits)
 *		bv_fill(bv, value)
 *              bv_set(bv, bitnum)
 *		bv_test(bv, bitnum)
 *		bv_nz(bv)
 *		bv_and(bv3, bv1, bv2)
 *		bv_or(bv3, bv1, bv2)
 *		bv_xor(bv3, bv1, bv2)
 *		bv_disp(title, bv)
 *
 *	\pre stdio.h
 *	\pre alloc.h
 *      \pre mem.h
 *
 *	ULTIMA MODIFICACION:
 *		LIBRARY.C, LIBRARY.H by A. J.-L. (1983)
 *		Traducción al castellano	 (1986)
 *		Extracción e implementación en PC-AT
 *		  (José Ramón Valverde Carrillo,  1989)
 *
 *	\author Andy Johnson-Laird.
 *		The programmer's CP/M handbook.
 *		Ed. Osborne / McGraw-Hill. 1983
 *
 *	COPYRIGHT:
 *		McGraw-Hill (1983)
 *
 *
 */

#include <stdio.h>

/** \struct _bv	
 *  Bit vector structure	
 */

struct _bv {
	unsigned bv_nobits;	/**< Number of bits in the vector     */
	char *bv_bits;		/**< Pointer to first bit (byte)	    */
	char *bv_end;		/**< Pointer to byte after bit vector */
};
typedef struct _bv *bit_vector;


/** bit_vector bv_new(bits)
 *      Create a new bit vector and set to zeros
 *  Reserve the memory space needed for the bit vector
 *  \param bits  size in bits of the vector
 *  \return a new bit vector, or NULL if there is not enough memory
 */
bit_vector bv_new(bits)
unsigned bits;			/* size in bits of the bit vector */
{
    extern void bv_fill();
    unsigned bytes;
    bit_vector bv;

    if ( !(bv = (bit_vector) malloc(sizeof(struct _bv))))
    	return NULL;
    bytes = (bits >> 3) + 1;
    if ( !(bv->bv_bits = (char *) malloc(bytes))
      return NULL;

    bv->bv_nobits = bits;			/* Fijación de longitud */
    bv->bv_end = bv->bv_bits + bytes;	/* Fijación de apuntador al final */

    bv_fill(bv, 0);				/* Rellenar vector con ceros */
    return bv;

} /* Final de bv_new */


/** void bv_fill(bv, value)  	
 * Flush the vector with a given value.
 *
 * This function sets the specified value in the given bit vector.
 * It only exists for consistency reasons, to isolate the module
 * from other standard functions like setmem()
 */
void bv_fill(bv, value)
bit_vector bv;			/* Pointer to a bit vector structure */
char value;			/* Value to set */
{
    /*     address      length     	    	   value */
    setmem(bv->bv_bits, (bv->bv_nobytes >> 3) + 1, value);
  }

/** void bv_set(bv, bitnum)		
 * set specified bit
 */
void bv_set(bv, bitnum)
bit_vector bv;			/**< \param Apuntador al vector de bits */
unsigned bitnum;		/**< \param Número de bit a fijar */

/* Results
   1 = OK
   0 = Error
*/
{
    unsigned byte_offset;	/* offset bytes in the vector */

    if ((byte_offset = bitnum >> 3) > bv -> bv_nobytes)
      return 0;		/* Attempt to overflow the vector */

    /*
     * Set the desired bit: we already know in which byte it is. The
     * bit number is easy to get using a logical AND with 0x07. Then
     * we just need to set it with a logical OR.
     */
     
    bv->bv_bits[byte_offset] |= (1 << (bitnum & 0x7));

    return 1;		/* OK */

} /* end of bv_set */

/**
 * bv_test(bv, bitnum)		
 * Test specified bit
 *
 * Esta función devuelve un valor que refleja el valor actual
 * del bit especificado. 
 */
bool bv_test(bv, bitnum)
/* Parámetros de entrada */
struct _bv *bv;			/* Apuntador al vector de bits */
unsigned bitnum;		/* Número de bit a examinar */

/* Parámetros de salida
   1 = OK!
   0 = Error
*/
{
    unsigned byte_offset;	/* Bytes de desplazamiento en el
				vector de bits */

    if ((byte_offset = bitnum >> 3) > bv -> bv_nobytes)
      return 0;		/* El número de bit sobrepasa el final del vector */

    /* Test de si está a uno el bit adecuado del vector. El byte de
    desplazamiento se ha calculado. El número del bit en el byte se
    calcula realizando la operación AND del número de bits con 0x07.
    El bit especificado se opera entonces (AND) con el vector */

    return bv -> bv_bits[byte_offset] & (1 << (bitnum & 0x7));

  } /* Final de bv_test */

/**
 * bv_nz(bv)		
 * Test de vector de bits distinto de cero 
 *
 * Esta función comprueba cada uno de los bytes
 * del vector especificado y vuelve indicando si existe algún
 * bit a uno en el vector.
 *
 */
int bv_nz(bv)
struct _bv *bv;		/* Apuntador al vector */

/* Parámetros de salida
   NZ = Uno de los bits del vector a uno
   0 = Todos los bits están a cero
*/

  {
    char *bits;		/* Apuntador a los bits del vector */

    bits = bv -> bv_bits;	/* Fijación de apuntador de trabajo */

    while (bits != bv -> bv_end)	/* Para vector de bits completo */
      {
	if (*bits++)		/* Si distinto de cero */
	  return *bits--;	/* Devolución del byte NZ */
      }
    return 0;		/* Indicador de vector cero */

  } /* Final de bv_nz */

/*===================================================================*/
void bv_and(bv3, bv1, bv2)			  /* bv3 = bv1 & bv2 */
/*===================================================================*/
/* Esta función realiza la operación AND entre los bytes de los vectores
   de bits 1 y 2 almacenando el resultado en el vector 3. */

/* Parámetros de entrada */
struct _bv *bv1;	/* Apuntador a vector de bits de entrada */
struct _bv *bv2;	/* Apuntador a vector de bits de entrada */

/* Parámetros de salida */
struct _bv *bv3;	/* Apuntador a vector de bits de salida */

  {
    char *bits1, *bits2, *bits3;	/* Apuntadores de trabajo a
					los vectores de bits */

    bits1 = bv1 -> bv_bits;		/* Inicialización de */
    bits2 = bv2 -> bv_bits;		/* apuntadores de trabajo */
    bits3 = bv3 -> bv_bits;

	/* Realización de la operación AND hasta el final de alguno
	de los vectores de bits */
    while ( (bits1 != bv1 -> bv_end) &&
	    (bits2 != bv2 -> bv_end) &&
	    (bits3 != bv3 -> bv_end) )
      {
	*bits3++ = *bits1++ & *bits2++; /* bv3 = bv1 & bv2 */
      }

  } /* Final de bv_and */

/*==========================================================*/
void bv_or(bv3, bv1, bv2)		/* bv3 = bv1 or bv2 */
/*==========================================================*/
/* Esta función realiza la operación OR entre los bytes de los vectores
   de bits 1 y 2 almacenando el resultado en el vector 3. */

/* Parámetros de entrada */
struct _bv *bv1;	/* Apuntador a vector de bits de entrada */
struct _bv *bv2;	/* Apuntador a vector de bits de entrada */

/* Parámetros de salida */
struct _bv *bv3;	/* Apuntador a vector de bits de salida */

  {
    char *bits1, *bits2, *bits3;	/* Apuntadores de trabajo a
					los vectores de bits */
    bits1 = bv1 -> bv_bits;		/* Inicialización de */
    bits2 = bv2 -> bv_bits;		/* apuntadores de trabajo */
    bits3 = bv3 -> bv_bits;

	/* Realización de la operación OR hasta el final de alguno
	de los vectores de bits */
    while ( (bits1 != bv1 -> bv_end) &&
	    (bits2 != bv2 -> bv_end) &&
	    (bits3 != bv3 -> bv_end) )
      {
	*bits3++ = *bits1++ | *bits2++; /* bv3 = bv1 or bv2 */
      }

  } /* Final de bv_or */

/*===================================================================*/
void bv_xor(bv3, bv1, bv2)		/* bv3 = bv1 xor bv2 */
/*===================================================================*/
/* Esta función realiza la operación XOR entre los bytes de los vectores
   de bits 1 y 2 almacenando el resultado en el vector 3. */

/* Parámetros de entrada */
struct _bv *bv1;	/* Apuntador a vector de bits de entrada */
struct _bv *bv2;	/* Apuntador a vector de bits de entrada */

/* Parámetros de salida */
struct _bv *bv3;	/* Apuntador a vector de bits de salida */

  {
    char *bits1, *bits2, *bits3;	/* Apuntadores de trabajo a
					los vectores de bits */

    bits1 = bv1 -> bv_bits;		/* Inicialización de apuntadores */
    bits2 = bv2 -> bv_bits;		/* de trabajo */
    bits3 = bv3 -> bv_bits;

	/* Realización de la operación XOR hasta el final de alguno
	de los vectores de bits */
    while ( (bits1 != bv1 -> bv_end) &&
	    (bits2 != bv2 -> bv_end) &&
	    (bits3 != bv3 -> bv_end) )
      {
	*bits3++ = *bits1++ ^ *bits2++; /* bv3 = bv1 xor bv2 */
      }

  } /* Final de bv_xor */

/*===================================================================*/
void bv_disp(title, bv)		/* Visualización de vectores de bits */
/*===================================================================*/
/* Esta función visualiza en hexadecimal el contenido del vector de bits
   especificado. Se utiliza normalmente para depuración de programas. */

/* Parámetros de entrada */
char *title;		/* Título para la visualización */
struct _bv *bv;

/* Parámetros de salida
   Ninguno.
*/

  {
    char *bits;			/* Apuntador de trabajo */
    unsigned byte_count;	/* Contador utilizado para formatear salida */
    unsigned bit_count;		/* Contador para procesamiento de
				los bits de un byte */
    char byte_value;		/* Valor a visualizar */

    fprintf(stderr, "\nBit Vector : %s", title);	/* Visualización
				de título */

    bits = bv -> bv_bits;	/* Fijación de apuntador de trabajo */
    byte_count = 0;		/* Inicialización del contador */

    while (bits != bv -> bv_end)	/* Para vector completo */
      {
	if ( (byte_count % 5) == 0)	/* Comprobación de nueva línea */
	  fprintf(stderr, "\n%4d : ", byte_count << 3);

        byte_value = *bits++;		/* Obtención del byte siguiente
					del vector */

	for (bit_count = 0; bit_count < 8; bit_count++)
	  {
	    /* Visualización del bit m s a la izquierda, a continuación
	    desplazamiento a la izquierda un bit */
	    if (bit_count == 4) fputc(' ', stderr); /* Separador */
	    fputc( (byte_value & 0x80) ? '1' : '0', stderr);
	    byte_value <<= 1;	/* Desplazamiento a la izquierda */
	  }
	fprintf(stderr, "    ");	/* Separador */

        byte_count++;		/* Actualización contador de bytes */
      }

  } /* Final de bv_disp */

/*	Final de P_BITS.C		*/


