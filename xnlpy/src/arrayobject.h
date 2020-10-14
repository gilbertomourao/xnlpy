/***************************************************************************
XNL: eXtended Numerical Library
2020, Gilberto Jose Guimaraes de Sousa Mourao

XNL is a library made by some students at UFMA (Universidade Federal do 
Maranhao) with the purpose of learning some numerical algorithms and 
share the acquired knowledge with students of the next generations.

                        About collaboration

Only UFMA students will be accepted as collaborators. This limitation was 
proposed by the creator of the library to prevent the focus from being 
diverted from learning.

If a student sees an error or feels that XNL has an incorrect algorithm, 
he should contact the library creator or the discipline teacher who will 
contact the library creator so that the student's complaint is studied. 
If there is really an error, the library will be updated and the student's 
name will appear in the acknowledgments (README.md).
***************************************************************************/

#ifndef XNLPY_H
	#error "arrayobject.h isn't intended to be included directly. Include xnl.h instead."
#endif

#ifndef ARRAYOBJECT_H
#define ARRAYOBJECT_H

#include <structmember.h>
#include <string.h>
#include <sys/types.h> /* for size_t */

/* Array type */
typedef struct xparray
{
	PyObject_HEAD /*declares ob_base of type PyObject*/
	double **data; /* pointer to the C array */
	int rows;
	int cols;
	struct xparray *aux; /* auxiliary object for 2D array assignment */
	int curr_row; /* current row, for a 2D array assignment */
	int curr_dim; /* current dimension, for 2D array slice */
	int row_step; /* row slice step for assignment */
} xparrayObject;

PyTypeObject xparrayType;

#define PyXParray_Check(obj) (strcmp(obj->ob_type->tp_name, "xnlpy.array"))

#endif /*ARRAYOBJECT_H*/