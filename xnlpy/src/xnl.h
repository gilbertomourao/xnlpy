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

#include <sys/time.h> /*struct timeval declaration*/
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <structmember.h>

/* Array type */
typedef struct 
{
	PyObject_HEAD /*declares ob_base of type PyObject*/
	double **data; /* pointer to the C array */
	int rows;
	int cols;
} xparrayObject;

PyTypeObject xparrayType;

/*Linear Algebra*/
xparrayObject *py_zeros(PyObject *, PyObject *);
xparrayObject *py_ones(PyObject *, PyObject *);
xparrayObject *py_eye(PyObject *, PyObject *);

/*Single Variable Calculus*/
PyObject *py_integral(PyObject *, PyObject *, PyObject *);