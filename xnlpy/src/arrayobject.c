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

#include "xnl.h"

#include <stdlib.h>

static PyMemberDef xparray_members[] = 
{
	/*{name, type, offset, flag, doc},*/
	{"rows", T_INT, offsetof(xparrayObject, rows), READONLY, "array rows"},
	{"cols", T_INT, offsetof(xparrayObject, cols), READONLY, "array columns"},
	{NULL} /*Sentinel*/
};

static int numParseArguments(double *array, Py_ssize_t size, PyObject *list)
{
	Py_ssize_t i;

	for (i = 0; i < size; i++)
	{
		PyObject *temp = PyList_GetItem(list, i); /*borrowed reference to returned value*/
		
		if (temp == NULL)
		{
			return 1;
		}

		/*Check if temp is numeric*/
		if (PyNumber_Check(temp) != 1)
		{
			return 1;
		}

		/*Convert temp to a python float*/
		temp = PyNumber_Float(temp);
		/*now to C double*/
		array[i] = PyFloat_AsDouble(temp);

		/*just to guarantee*/
		if (PyErr_Occurred())
			return 1;
	}

	return 0;
}

static double **ParseArguments(Py_ssize_t *ListSize, Py_ssize_t size, PyObject *args)
{
	static double **array;
	Py_ssize_t i;
	Py_ssize_t cols = 0;

	array = malloc(size * sizeof(double *));
	/*check the success of malloc call*/
	if (array == NULL) return NULL;

	for (i = 0; i < size; i++)
	{
		PyObject *temp = PyTuple_GetItem(args, i); /*borrowed reference to returned value*/
		/*checking if the argument is a list*/
		if (!PyList_CheckExact(temp))
		{
			PyErr_SetString(PyExc_TypeError, "Only lists are allowed as arguments.");
			return NULL;
		}
		/*It's a list, then continue*/
		/*Get the size of the list*/
		Py_ssize_t list_size = PyList_Size(temp);
		/*If there's more than one list, check if the other lists have all the same size of the first*/
		if (i > 0 && list_size != cols)
		{
			PyErr_SetString(PyExc_TypeError, "All the lists must have the same size.");
			return NULL;
		}
		cols = list_size;
		array[i] = malloc(list_size * sizeof(double));
		/*check the success of malloc call*/
		if (array[i] == NULL) return NULL;
		/*Now check if it's a list of numbers, not a list of lists or other objects*/
		if (numParseArguments(array[i], list_size, temp))
		{
			PyErr_SetString(PyExc_TypeError, "Only lists of numbers are allowed as arguments.");
			return NULL;
		}
	}

	*ListSize = cols;

	return array;
}

static int xparray_init(xparrayObject *self, PyObject *args, PyObject *kwds)
{
	Py_ssize_t TupleSize = PyTuple_Size(args);
	Py_ssize_t ListSize;

	if (!TupleSize)
	{
		PyErr_SetString(PyExc_TypeError, "You must supply at least one argument.");
		return -1;
	}

	self->data = ParseArguments(&ListSize, TupleSize, args);
	if (self->data == NULL)
		return -1;

	/*Passed*/
	self->rows = TupleSize;
	self->cols = ListSize;
	self->aux = NULL;
	self->ilow = 0;
	self->ihigh = 0;
	self->row_step = 0;
	self->tuplen = TupleSize;
	self->index_count = 0;

	return 0;
}

static PyObject *xparray_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	/*Will automatically call the init function*/
	xparrayObject *self;
	
	self = (xparrayObject *) type->tp_alloc(type, 0);

	if (self != NULL)
	{
		self->data = NULL;
		self->rows = 0;
		self->cols = 0;
		self->aux = NULL;
		self->ilow = 0;
		self->ihigh = 0;
		self->row_step = 0;
		self->tuplen = 0;
		self->index_count = 0;
	}

	return (PyObject *) self;
}

static void destroy(xparrayObject *self)
{
	int i;

	for (i = 0; i < self->rows; i++)
	{
		free(self->data[i]);
	}

	free(self->data);
}

static void xparray_dealloc(xparrayObject *self)
{
	destroy(self);
	Py_TYPE(self)->tp_free((PyObject *) self); /*to free the object's memory*/
}

/**
 * tp_as_method functions
 */

/**
 * print method
 *
 * I'm keeping this method here. It can print the xparray 
 * given a certain precision as parameter. It may be useful 
 * at some point(?).
 */
static PyObject *
xparray_print(xparrayObject *self, PyObject *args)
{
	int i = 0, j = 0, precision = 1;

	if (!PyArg_ParseTuple(args,"i",&precision))
	{
		return NULL;
	}

	if (precision < 0)
	{
		PyErr_SetString(PyExc_TypeError, "Precision must be a positive integer.");
		return NULL;
	}

	/*Print the array*/
	printf("array(");
	for (i = 0; i < self->rows; i++)
	{
		if (i) printf("%*s",6,"");
		printf("[");
		for (j = 0; j < self->cols; j++)
		{
			printf("%.*f, ", precision, self->data[i][j]);
		}
		printf("\b\b],");
		if (i < self->rows - 1)
			putchar('\n');
	}
	printf("\b)\n");

	return Py_BuildValue("");
}

static PyMethodDef xparray_methods[] = 
{
	{"print", (PyCFunction) xparray_print, METH_VARARGS, "Print the array"},

	{NULL} /*Sentinel*/
};

/**
 * tp_as_number functions
 */

xparrayObject *PyXParray_New(Py_ssize_t rows, Py_ssize_t cols)
{
	PyObject *argList = PyTuple_New(rows);
	Py_ssize_t i, j;

	for (i = 0; i < rows; i++)
	{
		PyObject *temp = PyList_New(cols);

		for (j = 0; j < cols; j++)
		{
			PyObject *zero = PyLong_FromLong(0);

			if (PyList_SetItem(temp, j, zero) == -1) /*steals the reference of zero*/
			{
				return NULL;
			}
		}

		PyTuple_SetItem(argList, i, temp); /*steals the reference of temp*/
	}

	/*call the array object*/
	xparrayObject *ret_array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (ret_array == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	return ret_array;
}

static PyObject *
xparray_add(PyObject *A, PyObject *B)
{
	double scalar;

	int A_scalar_flag = PyNumber_Check(A),
		B_scalar_flag = PyNumber_Check(B);

	if ((PyXParray_Check(A) && (A_scalar_flag != 1)) || (B_scalar_flag != 1 && PyXParray_Check(B)))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide two xnl arrays or one xnl array and a number.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array_A = NULL, *array_B = NULL;

	if (A_scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(A));
		array_B = (xparrayObject *) B;

		array_A = PyXParray_New(array_B->rows, array_B->cols);

		Py_ssize_t i, j;

		for (i = 0; i < array_B->rows; i++)
		{
			for (j = 0; j < array_B->cols; j++)
			{
				array_A->data[i][j] = scalar + array_B->data[i][j];
			}
		}

		return (PyObject *) array_A;
	}

	if (B_scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(B));
		array_A = (xparrayObject *) A;

		array_B = PyXParray_New(array_A->rows, array_A->cols);

		Py_ssize_t i, j;

		for (i = 0; i < array_A->rows; i++)
		{
			for (j = 0; j < array_A->cols; j++)
			{
				array_B->data[i][j] = array_A->data[i][j] + scalar;
			}
		}

		return (PyObject *) array_B;
	}

	array_A = (xparrayObject *) A;
	array_B = (xparrayObject *) B;

	/*Check if the element-wise operation can be done*/
	/**
	 * The two input xparrays must have compatible sizes (same dimensions or one of them = 1)
	 *
	 * Ex: 4x4 and 4x1 --> 4x4
	 *     1x2 and 3x5 --> incompatible
	 *     1x2 and 2x5 --> incompatible
	 *     2x1 and 2x5 --> 2x5
	 *
	 *  A->rows == B->rows or A->cols == B->cols when one dim of A and/or B = 1
	 *  dim(A) = dim(B)
	 */
	int a = array_A->rows == array_B->rows, 
		b = array_A->cols == array_B->cols,
		c = array_A->rows == 1,
		d = array_A->cols == 1,
		e = array_B->rows == 1,
		f = array_B->cols == 1;
	if ( !( (a && b) || (b && c) || (a && d) || (b && e) || (a && f) ) )
	{
		PyErr_SetString(PyExc_TypeError, "Incompatible arrays. Element-wise can be done "
										 "only if the two arrays are compatibles, i.e., "
										 "if both have the same dimensions or at least one "
										 "equivalent dimension when the other is one.");
		return NULL;
	}

	Py_ssize_t i, j;
	Py_ssize_t rows = (array_A->rows >= array_B->rows) ? array_A->rows : array_B->rows;
	Py_ssize_t cols = (array_A->cols >= array_B->cols) ? array_A->cols : array_B->cols;	

	xparrayObject *ret_array = PyXParray_New(rows, cols);

	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			ret_array->data[i][j] = array_A->data[i * !c][j * !d] + array_B->data[i * !e][j * !f];
		}
	}

	return (PyObject *) ret_array;
}

static PyObject *
xparray_sub(PyObject *A, PyObject *B)
{
	double scalar;

	int A_scalar_flag = PyNumber_Check(A),
		B_scalar_flag = PyNumber_Check(B);

	if ((PyXParray_Check(A) && (A_scalar_flag != 1)) || (B_scalar_flag != 1 && PyXParray_Check(B)))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide two xnl arrays or one xnl array and a number.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array_A = NULL, *array_B = NULL;

	if (A_scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(A));
		array_B = (xparrayObject *) B;

		array_A = PyXParray_New(array_B->rows, array_B->cols);

		Py_ssize_t i, j;

		for (i = 0; i < array_B->rows; i++)
		{
			for (j = 0; j < array_B->cols; j++)
			{
				array_A->data[i][j] = scalar - array_B->data[i][j];
			}
		}

		return (PyObject *) array_A;
	}

	if (B_scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(B));
		array_A = (xparrayObject *) A;

		array_B = PyXParray_New(array_A->rows, array_A->cols);

		Py_ssize_t i, j;

		for (i = 0; i < array_A->rows; i++)
		{
			for (j = 0; j < array_A->cols; j++)
			{
				array_B->data[i][j] = array_A->data[i][j] - scalar;
			}
		}

		return (PyObject *) array_B;
	}

	array_A = (xparrayObject *) A;
	array_B = (xparrayObject *) B;

	/*Check if the element-wise operation can be done*/
	/**
	 * The two input xparrays must have compatible sizes (same dimensions or one of them = 1)
	 *
	 * Ex: 4x4 and 4x1 --> 4x4
	 *     1x2 and 3x5 --> incompatible
	 *     1x2 and 2x5 --> incompatible
	 *     2x1 and 2x5 --> 2x5
	 *
	 *  A->rows == B->rows or A->cols == B->cols when one dim of A and/or B = 1
	 *  dim(A) = dim(B)
	 */
	int a = array_A->rows == array_B->rows, 
		b = array_A->cols == array_B->cols,
		c = array_A->rows == 1,
		d = array_A->cols == 1,
		e = array_B->rows == 1,
		f = array_B->cols == 1;
	if ( !( (a && b) || (b && c) || (a && d) || (b && e) || (a && f) ) )
	{
		PyErr_SetString(PyExc_TypeError, "Incompatible arrays. Element-wise can be done "
										 "only if the two arrays are compatibles, i.e., "
										 "if both have the same dimensions or at least one "
										 "equivalent dimension when the other is one.");
		return NULL;
	}

	Py_ssize_t i, j;
	Py_ssize_t rows = (array_A->rows >= array_B->rows) ? array_A->rows : array_B->rows;
	Py_ssize_t cols = (array_A->cols >= array_B->cols) ? array_A->cols : array_B->cols;	

	xparrayObject *ret_array = PyXParray_New(rows, cols);

	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			ret_array->data[i][j] = array_A->data[i * !c][j * !d] - array_B->data[i * !e][j * !f];
		}
	}

	return (PyObject *) ret_array;
}

static PyObject *
xparray_mul(PyObject *A, PyObject *B)
{
	double scalar;

	int A_scalar_flag = PyNumber_Check(A),
		B_scalar_flag = PyNumber_Check(B);

	if ((PyXParray_Check(A) && (A_scalar_flag != 1)) || (B_scalar_flag != 1 && PyXParray_Check(B)))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide two xnl arrays or one xnl array and a number.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array_A = NULL, *array_B = NULL;

	if (A_scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(A));
		array_B = (xparrayObject *) B;

		array_A = PyXParray_New(array_B->rows, array_B->cols);

		Py_ssize_t i, j;

		for (i = 0; i < array_B->rows; i++)
		{
			for (j = 0; j < array_B->cols; j++)
			{
				array_A->data[i][j] = scalar * array_B->data[i][j];
			}
		}

		return (PyObject *) array_A;
	}

	if (B_scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(B));
		array_A = (xparrayObject *) A;

		array_B = PyXParray_New(array_A->rows, array_A->cols);

		Py_ssize_t i, j;

		for (i = 0; i < array_A->rows; i++)
		{
			for (j = 0; j < array_A->cols; j++)
			{
				array_B->data[i][j] = array_A->data[i][j] * scalar;
			}
		}

		return (PyObject *) array_B;
	}

	array_A = (xparrayObject *) A;
	array_B = (xparrayObject *) B;

	/*Check if the element-wise operation can be done*/
	/**
	 * The two input xparrays must have compatible sizes (same dimensions or one of them = 1)
	 *
	 * Ex: 4x4 and 4x1 --> 4x4
	 *     1x2 and 3x5 --> incompatible
	 *     1x2 and 2x5 --> incompatible
	 *     2x1 and 2x5 --> 2x5
	 *
	 *  A->rows == B->rows or A->cols == B->cols when one dim of A and/or B = 1
	 *  dim(A) = dim(B)
	 */
	int a = array_A->rows == array_B->rows, 
		b = array_A->cols == array_B->cols,
		c = array_A->rows == 1,
		d = array_A->cols == 1,
		e = array_B->rows == 1,
		f = array_B->cols == 1;
	if ( !( (a && b) || (b && c) || (a && d) || (b && e) || (a && f) ) )
	{
		PyErr_SetString(PyExc_TypeError, "Incompatible arrays. Element-wise can be done "
										 "only if the two arrays are compatibles, i.e., "
										 "if both have the same dimensions or at least one "
										 "equivalent dimension when the other is one.");
		return NULL;
	}

	Py_ssize_t i, j;
	Py_ssize_t rows = (array_A->rows >= array_B->rows) ? array_A->rows : array_B->rows;
	Py_ssize_t cols = (array_A->cols >= array_B->cols) ? array_A->cols : array_B->cols;	

	xparrayObject *ret_array = PyXParray_New(rows, cols);

	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			ret_array->data[i][j] = array_A->data[i * !c][j * !d] * array_B->data[i * !e][j * !f];
		}
	}

	return (PyObject *) ret_array;
}

static PyObject *
xparray_matrix_mul(PyObject *A, PyObject *B)
{
	if (PyXParray_Check(A) || PyXParray_Check(B))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide two xparrays.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array_A = (xparrayObject *) A, 
				  *array_B = (xparrayObject *) B;

	/*Check if the multiplication can be done*/
	if (array_A->cols != array_B->rows)
	{
		PyErr_SetString(PyExc_TypeError, "The number of columns of A must be equal to the number of rows of B.");
		return NULL;
	}

	/*Passed, initialize a new array with zeros*/

	Py_ssize_t i, j;
	Py_ssize_t rows = array_A->rows;
	Py_ssize_t cols = array_B->cols;

	PyObject *argList = PyTuple_New(rows);

	for (i = 0; i < rows; i++)
	{
		PyObject *temp = PyList_New(cols);

		for (j = 0; j < cols; j++)
		{
			PyObject *zero = PyLong_FromLong(0);

			if (PyList_SetItem(temp, j, zero) == -1) /*steals the reference of zero*/
			{
				return NULL;
			}
		}

		PyTuple_SetItem(argList, i, temp); /*steals the reference of temp*/
	}

	/*call the array object*/
	xparrayObject *ret_array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (ret_array == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	/*Success on the creation, then do the multiplication*/
	Py_ssize_t k;

	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			for (k = 0; k < array_B->rows; k++)
			{
				ret_array->data[i][j] += array_A->data[i][k] * array_B->data[k][j];
			}
		}
	}

	/*Check if the output is a number.*/
	if (rows == 1 && cols == 1)
	{
		return Py_BuildValue("d", ret_array->data[0][0]);
	}

	return (PyObject *) ret_array;
}

static PyNumberMethods xparray_as_number = 
{
	(binaryfunc) xparray_add,          /* binaryfunc nb_add; */
    (binaryfunc) xparray_sub,          /* binaryfunc nb_subtract; */
    (binaryfunc) xparray_mul,          /* binaryfunc nb_multiply; */
    0,          		/* binaryfunc nb_remainder; */
    0,       			/* binaryfunc nb_divmod; */
    0,          		/* ternaryfunc nb_power; */
    0, 					/* unaryfunc nb_negative; */
    0, 					/* unaryfunc nb_positive; */
    0, 					/* unaryfunc nb_absolute; */
    0, 					/* inquiry nb_bool; */
    0,                  /* unaryfunc nb_invert; */
    0,                  /* binaryfunc nb_lshift; */
    0,                  /* binaryfunc nb_rshift; */
    0,                  /* binaryfunc nb_and; */
    0,                  /* binaryfunc nb_xor; */
    0,                  /* binaryfunc nb_or; */
    0,        			/* unaryfunc nb_int; */
    0,         			/* void *nb_reserved;  the slot formerly known as nb_long */
    0,        			/* unaryfunc nb_float; */
    0,                  /* binaryfunc nb_inplace_add; */
    0,                  /* binaryfunc nb_inplace_subtract; */
    0,                  /* binaryfunc nb_inplace_multiply; */
    0,                  /* binaryfunc nb_inplace_remainder; */
    0,                  /* ternaryfunc nb_inplace_power; */
    0,                  /* binaryfunc nb_inplace_lshift; */
    0,                  /* binaryfunc nb_inplace_rshift; */
    0,                  /* binaryfunc nb_inplace_and; */
    0,                  /* binaryfunc nb_inplace_xor; */
    0,                  /* binaryfunc nb_inplace_or; */
    0, 					/* binaryfunc nb_floor_divide; */
    0,					/* binaryfunc nb_true_divide; */
    0,					/* binaryfunc nb_inplace_floor_divide; */
    0,          		/* binaryfunc nb_inplace_true_divide; */
    0, 					/* unaryfunc nb_index; */
    (binaryfunc) xparray_matrix_mul, 					/* binaryfunc nb_matrix_multiply; */
    0,					/* binaryfunc nb_inplace_matrix_multiply; */
};

/**
 * 
 * tp_as_mapping functions
 * 
 */

/**
 * len() function
 *
 * Returns the largest dimension of the xparray. 
 *
 * e.g. 
 * >> A = xp.array([1,2,3],[4,5,6])
 * >> len(A)
 * 3
 * 
 */
static Py_ssize_t
xparray_length(xparrayObject *a)
{
	return (a ->rows >= a->cols) ? a->rows : a->cols;
}

/**
 * For index errors
 */
static PyObject *indexerr = NULL;

/**
 * valid_index()
 * 
 * To check if it is a valid index or not.
 * It follows the same implementation off listobject.c (see source code on github)
 */
static inline int 
valid_index(Py_ssize_t i, Py_ssize_t limit)
{
	return (size_t) i < (size_t) limit;
}

/**
 * xparray_item()
 *
 * Returns an object. There are the following cases:
 *
 * A[i] : returns a xparray with 1 row and cols columns if
 * 		  the original array A is a 2d array. If it is a 1d 
 * 		  array, then returns a float python object corresponding
 * 		  to the double in the position i.
 *
 * A[i][j] : returns a float pyton object corresponding to
 * 			 the double in the position [i][j].
 *
 * A[io:if:s][j] : Returns the column j from line io to line if, 
 * 				   with step s. If io = if, returns a float python 
 * 				   corresponding to the position [io][j].
 */
static PyObject *
xparray_item(xparrayObject *a, Py_ssize_t pos)
{
	Py_ssize_t limit = (a->rows == 1) ? a->cols : a->rows;
	
	if (!valid_index(pos, limit))
	{
		if (indexerr == NULL)
		{
			indexerr = PyUnicode_FromString("xparray index out of range");
			if (indexerr == NULL)
				return NULL;
		}
		PyErr_SetObject(PyExc_IndexError, indexerr);
		return NULL;
	}

	if (a->aux == NULL)
	{
		/*Reset parameters*/
		a->ilow = 0;
		a->ihigh = 0;
		a->row_step = 0;
		a->tuplen = 1;
	}

	if (a->rows == 1)
	{
		return Py_BuildValue("d", a->data[0][pos]); /*1D array*/
	}
	else
	{
		/*Returns a new array*/
		/*array[i]*/
		Py_ssize_t i, j;
		Py_ssize_t i_init = pos * (a->aux == NULL) + a->ilow * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
		Py_ssize_t i_end = pos * (a->aux == NULL) + a->ihigh * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
		Py_ssize_t j_init = 0 * (a->aux == NULL) + pos * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
		Py_ssize_t j_end = (a->cols - 1) * (a->aux == NULL) + pos * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
		Py_ssize_t tup_len = 0;

		PyObject *argList = PyTuple_New(a->tuplen);

		for (i = i_init; (a->row_step >= 0) ? i <= i_end : i >= i_end; i += (a->row_step == 0) ? 1 : a->row_step)
		{
			PyObject *temp = PyList_New((a->row_step == 0) ? a->cols : 1);
			Py_ssize_t j_list = 0;

			for (j = j_init; j <= j_end; j++)
			{
				xparrayObject *array = (a->aux == NULL) ? a : a->aux;
				PyObject *number = PyFloat_FromDouble(array->data[i][j]);

				if (PyList_SetItem(temp, j_list++, number) == -1) /*steals the reference of number*/
				{
					return NULL;
				}
			}

			PyTuple_SetItem(argList, tup_len++, temp); /*steals the reference of temp*/
		}

		/*call the array object*/
		xparrayObject *ret_array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

		Py_DECREF(argList);

		if (ret_array == NULL)
		{
			/*Will print the messages from xparray*/
			return NULL;
		}

		/*reset aux*/
		if (a->aux == NULL)
		{
			ret_array->aux = a;
			if (a->rows == 1) /*1D array*/
				ret_array->index_count = 2; /*forces index_count to 2*/
		} 
		else 
		{
			ret_array->aux = NULL;
			ret_array->index_count = 2;
		}

		ret_array->ilow = pos;
		ret_array->ihigh = pos;
		ret_array->tuplen = 1;
		/*keeps ret_array->row_step equals to the current value*/

		return (PyObject *) ret_array;
	}
}

/**
 * xparray_slice()
 *
 * Returns an object. There are the following cases:
 *
 * A[io:if:si] : returns a xparray with len rows and cols columns if
 * 		  		 the original array A is a 2d array, where len is the 
 * 		  		 number of elements on the interval [io,if], given a 
 * 		  		 step si. If it is a 1d array, then returns a xparray 
 * 		  		 with 1 row and len cols. If len = 1, then returns a 
 * 		  		 python float corresponding to the position [io].
 *
 * A[i][jo:jf:sj] : returns a xparray with 1 row and len columns, where
 * 					len is the number of elements on the interval [jo,jf], 
 * 					given a step sj. If len = 1, returns a python float 
 * 					corresponding to the position [i][jo].
 *
 * A[io:if:si][jo:jf:sj] : Returns a xparray with len_i rows and len_j cols, 
 * 						   where len_i is the number of elements on the 
 * 						   interval [io,if], given a step si, and len_j is 
 * 						   the number of elements on the interval [jo,jf], 
 * 						   given a step sj. If io=if and jo=jf, returns a 
 * 						   python float corresponding to the position [io][jo].
 */
static PyObject *
xparray_slice(xparrayObject *a, Py_ssize_t ilow, Py_ssize_t ihigh, Py_ssize_t step, Py_ssize_t limit)
{
	Py_ssize_t len = 0, iaux = ilow;

	do
	{
		iaux += step;
		len++;
	} while ((step > 0) ? iaux <= ihigh : iaux >= ihigh);

	if (len > limit)
	{
		PyErr_SetString(PyExc_IndexError, "xparray index out of range");
		return NULL;
	}

	if (a->aux == NULL)
	{
		a->ilow = ilow;
		a->ihigh = ihigh;
		a->row_step = step;
		a->tuplen = (a->rows == 1) ? 1: len;
	}
	else if ((a->ilow == a->ihigh) && (ilow == ihigh))
	{
		return Py_BuildValue("d", a->aux->data[a->ilow][ilow]);
	}

	Py_ssize_t i, j;
	Py_ssize_t i_init = ilow * (a->aux == NULL && a->rows > 1) + a->ilow * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
	Py_ssize_t i_end = ihigh * (a->aux == NULL && a->rows > 1) + a->ihigh * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
	Py_ssize_t j_init = 0 * (a->aux == NULL && a->rows > 1) + ilow * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
	Py_ssize_t j_end = (a->cols - 1) * (a->aux == NULL && a->rows > 1) + ihigh * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
	Py_ssize_t j_step = 0 * (a->aux == NULL && a->rows > 1) + step * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
	Py_ssize_t tup_len = 0;
	Py_ssize_t list_len;

	if (a->aux == NULL)
	{
		if (a->rows == 1)
			list_len = len;
		else 
			list_len = a->cols;
	} else 
		list_len = len;

	PyObject *argList = PyTuple_New(a->tuplen);

	for (i = i_init; (a->row_step >= 0) ? i <= i_end : i >= i_end; i += (a->row_step == 0) ? 1 : a->row_step)
	{
		PyObject *temp = PyList_New(list_len);
		Py_ssize_t j_list = 0;

		for (j = j_init; (j_step >= 0) ? j <= j_end : j >= j_end; j += (j_step == 0) ? 1 : j_step)
		{
			xparrayObject *array = (a->aux == NULL) ? a : a->aux;
			PyObject *number = PyFloat_FromDouble(array->data[i][j]);

			if (PyList_SetItem(temp, j_list++, number) == -1) /*steals the reference of number*/
			{
				return NULL;
			}
		}

		PyTuple_SetItem(argList, tup_len++, temp); /*steals the reference of temp*/
	}

	/*call the array object*/
	xparrayObject *ret_array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (ret_array == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	/*reset aux*/
	if (a->aux == NULL)
	{
		ret_array->aux = a;
		if (a->rows == 1) /*1D array*/
			ret_array->index_count = 2; /*forces index_count to 2*/
	} 
	else 
	{
		ret_array->aux = NULL;
		ret_array->index_count = 2;
	}

	ret_array->ilow = ilow;
	ret_array->ihigh = ihigh;
	ret_array->row_step = step;

	return (PyObject *) ret_array;

}

static PyObject *
xparray_subscript(xparrayObject* self, PyObject *item)
{
	if (self->index_count == 2)
	{
		PyErr_SetString(PyExc_IndexError, 
						"xparray doesn't work with more than one (1D array) or two (2D array) brackets.");
		return NULL;
	}

	if (PyIndex_Check(item)) /*returns 1 if item is an index integer (ex: [1]) and 0 otherwise*/
	{
		Py_ssize_t i;
		i = PyNumber_AsSsize_t(item, PyExc_IndexError); 
		/*returns item converted to a Py_ssize_t value it it can be interpreted as an integer.*/
		/*If fails, returns -1 and raises PyExc_IndexError*/
		if (i == -1 && PyErr_Occurred())
		{
			return NULL;
		}
		if (i < 0)
		{
			/* i += PyList_GET_SIZE(self) */ /*len(list)*/
			i += (self->rows == 1) ? self->cols : self->rows;
		}
		return xparray_item(self, i); 
	}
	else if (PySlice_Check(item))
	{
		Py_ssize_t start, stop, step;

		if (PySlice_Unpack(item, &start, &stop, &step) < 0)
		{
			return NULL;
		}

		if (step == 0)
		{
			PyErr_SetString(PyExc_IndexError, 
						"xparray slice index step can not be zero");
			return NULL;	
		}

		Py_ssize_t limit;

		if (self->aux == NULL)
		{
			limit = (self->rows == 1) ? self->cols : self->rows;
		}
		else 
		{
			limit = self->cols;
		}

		start += (start < 0) * limit;
		stop += (stop < 0) * limit;

		if (!valid_index(start, limit) || !valid_index(stop, limit))
		{
			PyErr_SetString(PyExc_IndexError, 
						"xparray index out of range");
			return NULL;
		}

		if (start > stop && step > 0)
		{
			PyErr_SetString(PyExc_IndexError, 
						"xparray start index can not be greater than the stop index when step is positive");
			return NULL;	
		}

		if (start < stop && step < 0)
		{
			PyErr_SetString(PyExc_IndexError, 
						"xparray stop index can not be greater than the start index when step is negative");
			return NULL;	
		}		

		return xparray_slice(self, start, stop, step, limit);
	}
	else 
	{
		PyErr_Format(PyExc_TypeError, 
					 "xparray indices must be integers or slices, not %.200s",
					 Py_TYPE(item)->tp_name);
		return NULL;
	}
}

/**
 * xparray_ass_item()
 *
 * Same of xparray_item, but with assignment.
 *
 * You can assign a xparray to a xparray or a number to a xparray. 
 * In the first case, be sure that both xparrays have the same dimensions. 
 * In the last case, the number will be assigned to all positions of the 
 * xparray argument.
 */
static int 
xparray_ass_item(xparrayObject *a, Py_ssize_t pos, PyObject *v)
{
	Py_ssize_t limit = (a->rows == 1) ? a->cols : a->rows;
	if (!valid_index(pos, limit))
	{
		PyErr_SetString(PyExc_IndexError, 
						"xparray assignment index out of range");
		return -1;
	}

	if (a->aux == NULL)
	{
		/*Reset parameters*/
		a->ilow = 0;
		a->ihigh = 0;
		a->row_step = 0;
		a->tuplen = 1;
	}

	/*Check if v is numeric*/
	if (PyNumber_Check(v) != 1)
	{
		if (!PyXParray_Check(v))
		{
			/*If the input is a xparray*/
			xparrayObject *array_v = (xparrayObject *) v;

			/*Check if the arrays' dimensions are equal*/
			if (a->tuplen == array_v->rows && a->cols == array_v->cols)
			{
				Py_ssize_t i, j;
				Py_ssize_t i_init = pos * (a->aux == NULL) + a->ilow * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
				Py_ssize_t i_end = pos * (a->aux == NULL) + a->ihigh * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
				Py_ssize_t j_init = 0 * (a->aux == NULL) + pos * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
				Py_ssize_t j_end = (a->cols - 1) * (a->aux == NULL) + pos * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
				Py_ssize_t i_v, j_v;

				for (i = i_init, i_v = 0; (a->row_step >= 0) ? i <= i_end : i >= i_end; i += (a->row_step == 0) ? 1 : a->row_step, i_v++)
				{
					for (j = j_init, j_v = 0; j <= j_end; j++, j_v++)
					{
						xparrayObject *array = (a->aux == NULL) ? a : a->aux;
						array->data[i][j] = array_v->data[i_v][j_v];
					}
				}

				return 0;
			}
			else 
			{
				PyErr_SetString(PyExc_ValueError, "When assign a xparray to another xparray, be sure that both have the same dimensions");
				return -1;
			}
		}

		PyErr_SetString(PyExc_TypeError, "You must assign a number to a xparray or a xparray to a xparray");

		return -1;
	}
	else
	{
		/*Convert v to a python float*/
		v = PyNumber_Float(v);
		/*now to C double*/
		double dv = PyFloat_AsDouble(v);
		
		Py_ssize_t i, j;
		Py_ssize_t i_init = pos * (a->aux == NULL && a->rows != 1) + a->ilow * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
		Py_ssize_t i_end = pos * (a->aux == NULL && a->rows != 1) + a->ihigh * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
		Py_ssize_t j_init = 0 * (a->aux == NULL && a->rows != 1) + pos * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
		Py_ssize_t j_end = (a->cols - 1) * (a->aux == NULL && a->rows != 1) + pos * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
		
		for (i = i_init; (a->row_step >= 0) ? i <= i_end : i >= i_end; i += (a->row_step == 0) ? 1 : a->row_step)
		{
			for (j = j_init; j <= j_end; j++)
			{
				xparrayObject *array = (a->aux == NULL) ? a : a->aux;
				array->data[i][j] = dv;
			}
		}

		/*just to guarantee*/
		if (PyErr_Occurred())
			return -1;
	}

	return 0;
}

/**
 * xparray_ass_slice()
 * 
 * Same of xparray_slice, but with assignment.
 *
 * You can assign a xparray to a xparray or a number to a xparray. 
 * In the first case, be sure that both xparrays have the same dimensions. 
 * In the last case, the number will be assigned to all positions of the 
 * xparray argument.
 */
static int 
xparray_ass_slice(xparrayObject *a, Py_ssize_t ilow, Py_ssize_t ihigh, Py_ssize_t step, Py_ssize_t limit, PyObject *v)
{
	Py_ssize_t len = 0, iaux = ilow;

	do
	{
		iaux += step;
		len++;
	} while ((step > 0) ? iaux <= ihigh : iaux >= ihigh);

	if (len > limit)
	{
		PyErr_SetString(PyExc_IndexError, "xparray index out of range");
		return -1;
	}

	if (a->aux == NULL)
	{
		a->ilow = ilow;
		a->ihigh = ihigh;
		a->row_step = step;
		a->tuplen = (a->rows == 1) ? 1: len;
	}

	/*Check if v is numeric*/
	if (PyNumber_Check(v) != 1)
	{
		if (!PyXParray_Check(v))
		{
			/*If the input is a xparray*/
			xparrayObject *array_v = (xparrayObject *) v;
			Py_ssize_t cols_len = (a->aux == NULL) ? a->cols : len;

			/*Check if the arrays' dimensions are equal*/
			if (a->tuplen == array_v->rows && cols_len == array_v->cols)
			{
				Py_ssize_t i, j;
				Py_ssize_t i_init = ilow * (a->aux == NULL) + a->ilow * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
				Py_ssize_t i_end = ihigh * (a->aux == NULL) + a->ihigh * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
				Py_ssize_t j_init = 0 * (a->aux == NULL) + ilow * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
				Py_ssize_t j_end = (a->cols - 1) * (a->aux == NULL) + ihigh * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
				Py_ssize_t j_step = 0 * (a->aux == NULL) + step * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
				Py_ssize_t i_v, j_v;

				for (i = i_init, i_v = 0; (a->row_step >= 0) ? i <= i_end : i >= i_end; i += (a->row_step == 0) ? 1 : a->row_step, i_v++)
				{
					for (j = j_init, j_v = 0; (j_step >= 0) ? j <= j_end : j >= j_end; j += (j_step == 0) ? 1 : j_step, j_v++)
					{
						xparrayObject *array = (a->aux == NULL) ? a : a->aux;
						array->data[i][j] = array_v->data[i_v][j_v];
					}
				}

				return 0;
			}
			else 
			{
				PyErr_SetString(PyExc_ValueError, "When assign a xparray to another xparray, be sure that both have the same dimensions");
				return -1;
			}
		}

		PyErr_SetString(PyExc_TypeError, "You must assign a number to a xparray or a xparray to a xparray");

		return -1;
	}
	else
	{
		/*Convert v to a python float*/
		v = PyNumber_Float(v);
		/*now to C double*/
		double dv = PyFloat_AsDouble(v);
		
		Py_ssize_t i, j;
		Py_ssize_t i_init = ilow * (a->aux == NULL) + a->ilow * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
		Py_ssize_t i_end = ihigh * (a->aux == NULL) + a->ihigh * (a->aux != NULL) + 0 * (a->aux == NULL && a->rows == 1);
		Py_ssize_t j_init = 0 * (a->aux == NULL) + ilow * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
		Py_ssize_t j_end = (a->cols - 1) * (a->aux == NULL) + ihigh * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
		Py_ssize_t j_step = 0 * (a->aux == NULL) + step * (a->aux != NULL || (a->aux == NULL && a->rows == 1));
	
		for (i = i_init; (a->row_step >= 0) ? i <= i_end : i >= i_end; i += (a->row_step == 0) ? 1 : a->row_step)
		{
			for (j = j_init; (j_step >= 0) ? j <= j_end : j >= j_end; j += (j_step == 0) ? 1 : j_step)
			{
				xparrayObject *array = (a->aux == NULL) ? a : a->aux;
				array->data[i][j] = dv;
			}
		}

		/*just to guarantee*/
		if (PyErr_Occurred())
			return -1;
	}

	return 0;
}

static int 
xparray_ass_subscript(xparrayObject *self, PyObject *item, PyObject *value)
{
	if (self->index_count == 2)
	{
		PyErr_SetString(PyExc_IndexError, 
						"xparray doesn't work with more than one (1D array) or two (2D array) brackets.");
		return -1;
	}

	if (PyIndex_Check(item))
	{
		Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		if (i == -1 && PyErr_Occurred())
			return -1;
		if (i < 0)
			i += (self->rows == 1) ? self->cols : self->rows;
		return xparray_ass_item(self, i, value);
	}
	else if (PySlice_Check(item))
	{
		Py_ssize_t start, stop, step;

		if (PySlice_Unpack(item, &start, &stop, &step) < 0)
		{
			return -1;
		}

		if (step == 0)
		{
			PyErr_SetString(PyExc_IndexError, 
						"xparray slice index step can not be zero");
			return -1;	
		}

		Py_ssize_t limit;

		if (self->aux == NULL)
		{
			limit = (self->rows == 1) ? self->cols : self->rows;
		}
		else 
		{
			limit = self->cols;
		}

		start += (start < 0) * limit;
		stop += (stop < 0) * limit;

		if (!valid_index(start, limit) || !valid_index(stop, limit))
		{
			PyErr_SetString(PyExc_IndexError, 
						"xparray index out of range");
			return -1;
		}

		if (start > stop && step > 0)
		{
			PyErr_SetString(PyExc_IndexError, 
						"xparray start index can not be greater than the stop index when step is positive");
			return -1;	
		}

		if (start < stop && step < 0)
		{
			PyErr_SetString(PyExc_IndexError, 
						"xparray stop index can not be greater than the start index when step is negative");
			return -1;	
		}		

		return xparray_ass_slice(self, start, stop, step, limit, value);
	}
	else
	{
		PyErr_Format(PyExc_TypeError, 
					 "xparray indices must be integers or slices, not %.200s",
					 Py_TYPE(item)->tp_name);
		return -1;
	}
}

static PyMappingMethods xparray_as_mapping =
{
	(lenfunc) xparray_length,
	(binaryfunc) xparray_subscript,
	(objobjargproc) xparray_ass_subscript
};

/**
 * repr and str functions
 */

static PyObject *
xparray_repr(xparrayObject *self)
{
	Py_ssize_t i;
	PyObject *s; /*object representation of double*/
	_PyUnicodeWriter writer;

	i = Py_ReprEnter((PyObject *) self);
	if (i != 0)
	{
		return i > 0 ? PyUnicode_FromString("xparray([...])") : NULL;
	}

	_PyUnicodeWriter_Init(&writer);
	writer.overallocate = 1;
	/**
	 * First lines: 
	 * "xparray([" + "1" + ", 2" * (cols - 1) + "],\n"
	 * "        [" + "1" + ", 2" * (cols - 1) + "],\n"
	 * ...
	 * Last line:
	 * "        [" + "1" + ", 2" * (cols - 1) + "])"
	 */
	writer.min_length = self->rows * (9 + 1 + (2 + 1) * (self->cols - 1) + 1 + 1 + 1) - 1;

	if (_PyUnicodeWriter_WriteASCIIString(&writer, "xparray([", 9) < 0) /*returns -1 on error*/
		goto error; /*original implementation*/

	/*Do repr() on each element for rows and cols*/
	int i_row, j_col;

	for (i_row = 0; i_row < self->rows; i_row++)
	{
		for (j_col = 0; j_col < self->cols; j_col++)
		{
			if (j_col > 0)
			{
				if (_PyUnicodeWriter_WriteASCIIString(&writer, ", ", 2) < 0)
					goto error;
			}

			s = PyObject_Repr(Py_BuildValue("d", self->data[i_row][j_col]));
			if (s == NULL)
				goto error;

			if(_PyUnicodeWriter_WriteStr(&writer, s) < 0)
			{
				Py_DECREF(s);
				goto error;
			}

			Py_DECREF(s);
		}

		if (self->rows > 1 && i_row < self->rows - 1)
		{
			if (_PyUnicodeWriter_WriteASCIIString(&writer, "],\n        [", 12) < 0)
				goto error;
		}
		else
		{
			/*Last line*/
			if (_PyUnicodeWriter_WriteASCIIString(&writer, "])", 2) < 0)
				goto error;
		}	
	}

	writer.overallocate = 0;

	Py_ReprLeave((PyObject *) self);
	return _PyUnicodeWriter_Finish(&writer);

error: 
	_PyUnicodeWriter_Dealloc(&writer);
	Py_ReprLeave((PyObject *) self);
	return NULL;
}

/**
 * 
 */

PyTypeObject xparrayType = 
{
	PyVarObject_HEAD_INIT(NULL, 0) /*initialize ob_base*/
	.tp_name = "xnlpy.array", /*module.name is compatible with pydoc and pickle*/
	.tp_doc = "Array object used in xnlpy",
	.tp_basicsize = sizeof(xparrayObject),
	.tp_itemsize = 0, /*non zero only when dealing with objects with variable size*/
	.tp_flags = Py_TPFLAGS_DEFAULT, 
	.tp_new = xparray_new, /*enable object creation*/
	.tp_init = (initproc) xparray_init,
	.tp_new = xparray_new,
	.tp_dealloc = (destructor) xparray_dealloc,
	.tp_members = xparray_members,
	.tp_methods = xparray_methods,
	.tp_as_number = &xparray_as_number,
	.tp_as_mapping = &xparray_as_mapping,
	.tp_repr = (reprfunc) xparray_repr
};