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
	self->curr_row = 0;
	self->curr_dim = 0;

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
		self->curr_row = 0;
		self->curr_dim = 0;
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
/*
static PyObject *xparray_read(xparrayObject *self, PyObject *args)
{
	int i = 0, j = 0;
	double value;
	PyObject *retval;

	if (!PyArg_ParseTuple(args, "ii", &i, &j))
	{
		return NULL;
	}

	if (i < 0 || j < 0 || i > self->rows - 1 || j > self->cols - 1)
	{
		PyErr_SetString(PyExc_TypeError, "Invalid indices.");
		return NULL;
	}

	value = self->data[i][j];
	retval = Py_BuildValue("d", value);

	return retval;
}

static PyObject *xparray_write(xparrayObject *self, PyObject *args)
{
	int i, j;
	double value;

	if (!PyArg_ParseTuple(args, "iid", &i, &j, &value))
	{
		return NULL;
	}

	if (i < 0 || j < 0 || i > self->rows - 1 || j > self->cols - 1)
	{
		PyErr_SetString(PyExc_TypeError, "Invalid indices.");
		return NULL;
	}

	self->data[i][j] = value;

	return Py_BuildValue("");
}
*/
static PyObject *xparray_print(xparrayObject *self, PyObject *args)
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
	/*
	{"read", (PyCFunction) xparray_read, METH_VARARGS, "Return the element at (row, column)"},
	{"write", (PyCFunction) xparray_write, METH_VARARGS, "Set the element at (row, column)"},
	*/
	{"print", (PyCFunction) xparray_print, METH_VARARGS, "Print the array"},

	{NULL} /*Sentinel*/
};

/**
 * tp_as_number functions
 */

static PyObject *xparray_mul(PyObject *A, PyObject *B)
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

		/*Now create the array ([scalar])*/
		PyObject *argList = PyTuple_New(1);
		PyObject *temp = PyList_New(1);
		PyObject *py_scalar = PyFloat_FromDouble(scalar);

		if (PyList_SetItem(temp, 0, py_scalar) == -1) /*steals the reference of py_scalar*/
		{
			return NULL;
		}

		PyTuple_SetItem(argList, 0, temp); /*steals the reference of temp*/

		array_A = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);
		array_B = (xparrayObject *) B;
	}
	else if (B_scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(B));

		/*Now create the array ([scalar])*/
		PyObject *argList = PyTuple_New(1);
		PyObject *temp = PyList_New(1);
		PyObject *py_scalar = PyFloat_FromDouble(scalar);

		if (PyList_SetItem(temp, 0, py_scalar) == -1) /*steals the reference of py_scalar*/
		{
			return NULL;
		}

		PyTuple_SetItem(argList, 0, temp); /*steals the reference of temp*/

		array_B = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);
		array_A = (xparrayObject *) A;	
	}
	else /*arrays multiplication*/
	{
		array_A = (xparrayObject *) A;
		array_B = (xparrayObject *) B;

		/*Check if the multiplication can be done*/
		if (array_A->cols != array_B->rows)
		{
			PyErr_SetString(PyExc_TypeError, "The number of columns of A must be equal to the number of rows of B.");
			return NULL;
		}
	}

	/*Passed, initialize a new array with zeros*/

	int i, j;
	int rows = (A_scalar_flag == 1) ? array_B->rows : array_A->rows;
	int cols = (B_scalar_flag == 1) ? array_A->cols : array_B->cols;

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
	int k;

	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			if (!A_scalar_flag && !B_scalar_flag)
			{
				for (k = 0; k < array_B->rows; k++)
				{
					ret_array->data[i][j] += array_A->data[i][k] * array_B->data[k][j];
				}
			}
			else
			{
				ret_array->data[i][j] += (A_scalar_flag == 1) ? array_A->data[0][0] * array_B->data[i][j] : array_A->data[i][j] * array_B->data[0][0];
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

static PyObject *xparray_add(PyObject *A, PyObject *B)
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

		/*Now create the array ([scalar])*/
		PyObject *argList = PyTuple_New(1);
		PyObject *temp = PyList_New(1);
		PyObject *py_scalar = PyFloat_FromDouble(scalar);

		if (PyList_SetItem(temp, 0, py_scalar) == -1) /*steals the reference of py_scalar*/
		{
			return NULL;
		}

		PyTuple_SetItem(argList, 0, temp); /*steals the reference of temp*/

		array_A = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);
		array_B = (xparrayObject *) B;
	}
	else if (B_scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(B));

		/*Now create the array ([scalar])*/
		PyObject *argList = PyTuple_New(1);
		PyObject *temp = PyList_New(1);
		PyObject *py_scalar = PyFloat_FromDouble(scalar);

		if (PyList_SetItem(temp, 0, py_scalar) == -1) /*steals the reference of py_scalar*/
		{
			return NULL;
		}

		PyTuple_SetItem(argList, 0, temp); /*steals the reference of temp*/

		array_B = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);
		array_A = (xparrayObject *) A;	
	}
	else /*arrays multiplication*/
	{
		array_A = (xparrayObject *) A;
		array_B = (xparrayObject *) B;

		/*Check if the multiplication can be done*/
		if (array_A->rows != array_B->rows || array_A->cols != array_B->cols)
		{
			PyErr_SetString(PyExc_TypeError, "The parameters A and B must have the same dimensions to perform addition.\n");
			return NULL;
		}
	}

	/*Passed, initialize a new array with zeros*/

	int i, j;
	int rows = (A_scalar_flag == 1) ? array_B->rows : array_A->rows;
	int cols = (B_scalar_flag == 1) ? array_A->cols : array_B->cols;

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

	/*Success on the creation, then do the addition*/
	for (i = 0; i < ret_array->rows; i++)
	{
		for (j = 0; j < ret_array->cols; j++)
		{
			ret_array->data[i][j] = array_A->data[i*(!A_scalar_flag)][j*(!A_scalar_flag)] + array_B->data[i*(!B_scalar_flag)][j*(!B_scalar_flag)];
		}
	}

	return (PyObject *) ret_array;
}

static PyObject *xparray_sub(PyObject *A, PyObject *B)
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

		/*Now create the array ([scalar])*/
		PyObject *argList = PyTuple_New(1);
		PyObject *temp = PyList_New(1);
		PyObject *py_scalar = PyFloat_FromDouble(scalar);

		if (PyList_SetItem(temp, 0, py_scalar) == -1) /*steals the reference of py_scalar*/
		{
			return NULL;
		}

		PyTuple_SetItem(argList, 0, temp); /*steals the reference of temp*/

		array_A = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);
		array_B = (xparrayObject *) B;
	}
	else if (B_scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(B));

		/*Now create the array ([scalar])*/
		PyObject *argList = PyTuple_New(1);
		PyObject *temp = PyList_New(1);
		PyObject *py_scalar = PyFloat_FromDouble(scalar);

		if (PyList_SetItem(temp, 0, py_scalar) == -1) /*steals the reference of py_scalar*/
		{
			return NULL;
		}

		PyTuple_SetItem(argList, 0, temp); /*steals the reference of temp*/

		array_B = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);
		array_A = (xparrayObject *) A;	
	}
	else /*arrays multiplication*/
	{
		array_A = (xparrayObject *) A;
		array_B = (xparrayObject *) B;

		/*Check if the multiplication can be done*/
		if (array_A->rows != array_B->rows || array_A->cols != array_B->cols)
		{
			PyErr_SetString(PyExc_TypeError, "The parameters A and B must have the same dimensions to perform addition.\n");
			return NULL;
		}
	}

	/*Passed, initialize a new array with zeros*/

	int i, j;
	int rows = (A_scalar_flag == 1) ? array_B->rows : array_A->rows;
	int cols = (B_scalar_flag == 1) ? array_A->cols : array_B->cols;

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

	/*Success on the creation, then do the subtraction*/
	for (i = 0; i < ret_array->rows; i++)
	{
		for (j = 0; j < ret_array->cols; j++)
		{
			ret_array->data[i][j] = array_A->data[i*(!A_scalar_flag)][j*(!A_scalar_flag)] - array_B->data[i*(!B_scalar_flag)][j*(!B_scalar_flag)];
		}
	}

	return (PyObject *) ret_array;
}

static PyNumberMethods xparray_as_number = 
{
	xparray_add,          /*nb_add*/
    xparray_sub,          /*nb_subtract*/
    xparray_mul,          /*nb_multiply*/
    0, 					/*nb_divide*/
    0,          		/*nb_remainder*/
    0,       			/*nb_divmod*/
    0,          		/*nb_power*/
    0, 					/*nb_negative*/
    0, 					/*nb_positive*/
    0, 					/*nb_absolute*/
    0, 					/*nb_nonzero*/
    0,                  /*nb_invert*/
    0,                  /*nb_lshift*/
    0,                  /*nb_rshift*/
    0,                  /*nb_and*/
    0,                  /*nb_xor*/
    0,                  /*nb_or*/
    0,       			/*nb_coerce*/
    0,        			/*nb_int*/
    0,         			/*nb_long*/
    0,        			/*nb_float*/
    0,                  /* nb_oct */
    0,                  /* nb_hex */
    0,                  /* nb_inplace_add */
    0,                  /* nb_inplace_subtract */
    0,                  /* nb_inplace_multiply */
    0,                  /* nb_inplace_divide */
    0,                  /* nb_inplace_remainder */
    0,                  /* nb_inplace_power */
    0,                  /* nb_inplace_lshift */
    0,                  /* nb_inplace_rshift */
    0,                  /* nb_inplace_and */
    0,                  /* nb_inplace_xor */
    0,                  /* nb_inplace_or */
    0, 					/* nb_floor_divide */
    0,          		/* nb_true_divide */
};

/**
 * tp_as_mapping functions
 */

static Py_ssize_t
xparray_length(xparrayObject *a)
{
	return (a ->rows >= a->cols) ? a->rows : a->cols;
}

static PyObject *indexerr = NULL;

static inline int 
valid_index(Py_ssize_t i, Py_ssize_t limit)
{
	return (size_t) i < (size_t) limit;
}

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
	if (a->rows == 1)
	{
		return Py_BuildValue("d", a->data[0][pos]); /*1D array*/
	}
	else
	{
		/*Returns a new array*/
		/*array[i]*/
		int i, j;
		int rows = 1;
		int cols = a->cols;

		PyObject *argList = PyTuple_New(rows);

		for (i = 0; i < rows; i++)
		{
			PyObject *temp = PyList_New(cols);

			for (j = 0; j < cols; j++)
			{
				PyObject *number = PyFloat_FromDouble(a->data[pos][j]);

				if (PyList_SetItem(temp, j, number) == -1) /*steals the reference of number*/
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

		ret_array->curr_row = pos;
		ret_array->aux = a;
		ret_array->curr_dim = 1;

		return (PyObject *) ret_array;
	}
}

static PyObject *
xparray_slice(xparrayObject *a, Py_ssize_t ilow, Py_ssize_t ihigh, Py_ssize_t step)
{
	Py_ssize_t len = 0, iaux = ilow;

	do
	{
		iaux += step;
		len++;
	} while ((step > 0) ? iaux <= ihigh : iaux >= ihigh);

	if (a->curr_dim == 1)
	{
		if (len > a->cols)
		{
			PyErr_SetString(PyExc_IndexError, "xparray index out of range");
			return NULL;
		}

		if (len == 1)
		{
			return Py_BuildValue("d", a->data[0][ilow]);
		}

		/*1D array*/
		int i, j, list_len = 0;
		int rows = a->curr_row;
		int cols = len;

		PyObject *argList = PyTuple_New(rows);
		
		for (i = 0; i < rows; i++)
		{
			PyObject *temp = PyList_New(cols);

			for (j = ilow; (step > 0) ? j <= ihigh : j >= ihigh; j+=step)
			{
				PyObject *number = PyFloat_FromDouble(a->data[i][j]);

				if (PyList_SetItem(temp, list_len++, number) == -1) /*steals the reference of number*/
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

		return (PyObject *) ret_array;		
	}
	else
	{
		if (len > a->rows)
		{
			PyErr_SetString(PyExc_IndexError, "xparray index out of range");
			return NULL;
		}

		/*2D array*/
		Py_ssize_t i, j, tup_len = 0;
		Py_ssize_t rows = len;
		Py_ssize_t cols = a->cols;

		PyObject *argList = PyTuple_New(rows);

		for (i = ilow; (step > 0) ? i <= ihigh : i >= ihigh; i+=step)
		{
			PyObject *temp = PyList_New(cols);

			for (j = 0; j < cols; j++)
			{
				PyObject *number = PyFloat_FromDouble(a->data[i][j]);

				if (PyList_SetItem(temp, j, number) == -1) /*steals the reference of number*/
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

		ret_array->aux = a;
		ret_array->curr_row = ilow;
		ret_array->curr_dim = 1;
		ret_array->row_step = step;

		return (PyObject *) ret_array;
	}

}

static PyObject *
xparray_subscript(xparrayObject* self, PyObject *item)
{
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

/*
		Py_ssize_t length = (self->curr_dim == 1) ? self->cols : self->rows;

		slicelength = PySlice_AdjustIndices(length, &start, &stop, step);
*/	
		/*Now start and stop will be adjusted to the max/min values*/

		Py_ssize_t limit = (self->curr_dim == 1) ? self->cols : self->rows;

		start += (start < 0) * ( (self->rows == 1) ? self->cols : self->rows );
		stop += (stop < 0) * ( (self->rows == 1) ? self->cols : self->rows );

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

		return xparray_slice(self, start, stop, step);
	}
	else 
	{
		PyErr_Format(PyExc_TypeError, 
					 "xparray indices must be integers or slices, not %.200s",
					 Py_TYPE(item)->tp_name);
		return NULL;
	}
}

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

	/*Check if v is numeric*/
	if (PyNumber_Check(v) != 1)
	{
		if (!PyXParray_Check(v) && a->rows > 1)
		{
			/*If the input is a xparray*/
			xparrayObject *array_v = (xparrayObject *) v;

			/*Check if the arrays' dimensions are equal*/
			if (array_v->rows == 1 && a->cols == array_v->cols)
			{
				int i,j;

				for (i = 0; i < array_v->rows; i++)
				{
					for (j = 0; j < array_v->cols; j++)
					{
						a->data[pos][j] = array_v->data[i][j];
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

		PyErr_SetString(PyExc_TypeError, "You must assign a number to a xparray element or a xparray to a xparray");

		return -1;
	}
	else
	{
		/*Convert v to a python float*/
		v = PyNumber_Float(v);
		/*now to C double*/
		if (a->rows == 1)
		{
			a->aux->data[a->curr_row][pos] = PyFloat_AsDouble(v);
		}
		else 
		{
			PyErr_SetString(PyExc_TypeError, "cannot assign a number to a xparray");
			return -1;
		}

		/*just to guarantee*/
		if (PyErr_Occurred())
			return -1;
	}

	return 0;
}

static int 
xparray_ass_slice(xparrayObject *a, Py_ssize_t ilow, Py_ssize_t ihigh, Py_ssize_t step, PyObject *v)
{
	Py_ssize_t len = 0, iaux = ilow;

	do
	{
		iaux += step;
		len++;
	} while ((step > 0) ? iaux <= ihigh : iaux >= ihigh);

	if (a->curr_dim == 1)
	{
		if (len > a->cols)
		{
			PyErr_SetString(PyExc_IndexError, "xparray index out of range");
			return -1;
		}

		if (PyNumber_Check(v) != 1)
		{
			if (!PyXParray_Check(v))
			{
				/*If the input is a xparray*/
				xparrayObject *array_v = (xparrayObject *) v;

				/*Check if the arrays' dimensions are equal*/
				if (a->rows == array_v->rows && len == array_v->cols)
				{
					int i, j, k;

					for (i = 0; i < a->rows; i++)
					{
						for (j = ilow, k = 0; (step > 0) ? j <= ihigh : j >= ihigh; j += step, k++)
						{
							int row_index = a->curr_row + i*a->row_step;
							row_index += (row_index < 0) ? a->aux->rows: 0;
							a->aux->data[row_index][j] = array_v->data[i][k];
						}
					}
				}
				else 
				{
					PyErr_SetString(PyExc_ValueError, "When assign a xparray to another xparray, be sure that both have the same dimensions");
					return -1;
				}

				return 0;
			}

			PyErr_SetString(PyExc_TypeError, "You must assign a number to a xparray element or a xparray to a xparray");

			return -1;
		}

		/*v is a number, then assign it to all xparray elements*/
		/*Convert v to a python float*/
		v = PyNumber_Float(v);
		double val_v = PyFloat_AsDouble(v);

		int i, j;

		for (i = 0; i < a->rows; i++)
		{
			for (j = ilow; (step > 0) ? j <= ihigh : j >= ihigh; j += step)
			{
				int row_index = a->curr_row + i*a->row_step;
				row_index += (row_index < 0) ? a->aux->rows: 0;
				a->aux->data[row_index][j] = val_v;
			}
		}
	}
	else 
	{
		if (len > a->rows)
		{
			PyErr_SetString(PyExc_IndexError, "xparray index out of range");
			return -1;
		}

		if (PyNumber_Check(v) != 1)
		{
			if (!PyXParray_Check(v))
			{
				/*If the input is a xparray*/
				xparrayObject *array_v = (xparrayObject *) v;

				/*Check if the arrays' dimensions are equal*/
				if (len == array_v->rows && a->cols == array_v->cols)
				{
					int i, j, k;

					for (i = ilow, k = 0; (step > 0) ? i <= ihigh : i >= ihigh; i+=step, k++)
					{
						for (j = 0; j < a->cols;j++)
						{
							a->data[i][j] = array_v->data[k][j];
						}
					}
				}
				else 
				{
					PyErr_SetString(PyExc_ValueError, "When assign a xparray to another xparray, be sure that both have the same dimensions");
					return -1;
				}

				return 0;
			}

			PyErr_SetString(PyExc_TypeError, "You must assign a number to a xparray element or a xparray to a xparray");

			return -1;	
		}

		/*v is a number, then assign it to all xparray elements*/
		/*Convert v to a python float*/
		v = PyNumber_Float(v);
		double val_v = PyFloat_AsDouble(v);

		int i, j;

		for (i = ilow; (step > 0) ? i <= ihigh : i >= ihigh; i+=step)
		{
			for (j = 0; j < a->cols;j++)
			{
				a->data[i][j] = val_v;
			}
		}			
	}

	return 0;
}

static int 
xparray_ass_subscript(xparrayObject *self, PyObject *item, PyObject *value)
{
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

		Py_ssize_t limit = (self->curr_dim == 1) ? self->cols : self->rows;

		start += (start < 0) * ( (self->rows == 1) ? self->cols : self->rows );
		stop += (stop < 0) * ( (self->rows == 1) ? self->cols : self->rows );

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

		return xparray_ass_slice(self, start, stop, step, value);
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
	.tp_as_mapping = &xparray_as_mapping
};