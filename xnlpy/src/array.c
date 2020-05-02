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

	array = malloc(size * sizeof(double));
	/*check the success of malloc call*/
	if (array == NULL) return NULL;

	for (i = 0; i < size; i++)
	{
		PyObject *temp = PyTuple_GetItem(args, i); /*borrowed reference to returned value*/
		/*checking if the argument is a list*/
		if (!PyList_CheckExact(temp))
		{
			PyErr_SetString(PyExc_TypeError, "Only list are allowed as arguments.");
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
	Py_ssize_t *ListSize = malloc(sizeof(ListSize));

	if (!TupleSize)
	{
		PyErr_SetString(PyExc_TypeError, "You must supply at least one argument.");
		return -1;
	}

	self->data = ParseArguments(ListSize, TupleSize, args);
	if (self->data == NULL)
		return -1;

	/*Passed*/
	self->rows = TupleSize;
	self->cols = *ListSize;

	free(ListSize);

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

static PyMemberDef xparray_members[] = 
{
	/*{name, type, offset, flag, doc},*/
	{"rows", T_INT, offsetof(xparrayObject, rows), READONLY, "array rows"},
	{"cols", T_INT, offsetof(xparrayObject, cols), READONLY, "array columns"},
	{NULL} /*Sentinel*/
};

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
	int success = 1;
	PyObject *retval;

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
	retval = Py_BuildValue("i", success);

	return retval;
}

static PyObject *xparray_print(xparrayObject *self, PyObject *args)
{
	int i = 0, j = 0, precision = 1;
	int success = 1;
	PyObject *retval;

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

	retval = Py_BuildValue("i", success);

	return retval;
}

xparrayObject *py_dot(xparrayObject *self, PyObject *args)
{
	PyObject *array_arg = NULL;
	double scalar;

	if (!PyArg_ParseTuple(args, "O", &array_arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(array_arg);

	if (PyXParray_Check(array_arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "The argument's type must be xnlpy.array or numeric.");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		array_arg = PyNumber_Float(array_arg);
		scalar = PyFloat_AsDouble(array_arg);

		/*Now create the array ([scalar])*/
		PyObject *argList = PyTuple_New(1);
		PyObject *temp = PyList_New(1);
		PyObject *py_scalar = PyFloat_FromDouble(scalar);

		if (PyList_SetItem(temp, 0, py_scalar) == -1) /*steals the reference of py_scalar*/
		{
			return NULL;
		}

		PyTuple_SetItem(argList, 0, temp); /*steals the reference of temp*/

		array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);
	}
	else 
	{
		array = (xparrayObject *)array_arg;

		/*Check if the multiplication can be done*/
		if (self->cols != array->rows)
		{
			PyErr_SetString(PyExc_TypeError, "The number of columns of A must be equal to the number of rows of B.");
			return NULL;
		}
	}

	/*Passed, initialize a new array with zeros*/
	PyObject *argList = PyTuple_New(self->rows);

	int i, j;
	int cols = (scalar_flag == 1) ? self->cols : array->cols;

	for (i = 0; i < self->rows; i++)
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

	for (i = 0; i < self->rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			for (k = 0; k < array->rows; k++)
			{
				ret_array->data[i][j] += self->data[i][k + j*scalar_flag] * array->data[k][j*(1 - scalar_flag)];
			}
		}
	}

	return ret_array;
}

xparrayObject *py_plus(xparrayObject *self, PyObject *args)
{
	PyObject *array_arg = NULL;
	double scalar;

	if (!PyArg_ParseTuple(args, "O", &array_arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(array_arg);

	if (PyXParray_Check(array_arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "The argument's type must be xnlpy.array or numeric.");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		array_arg = PyNumber_Float(array_arg);
		scalar = PyFloat_AsDouble(array_arg);

		/*Now create the array ([scalar])*/
		PyObject *argList = PyTuple_New(1);
		PyObject *temp = PyList_New(1);
		PyObject *py_scalar = PyFloat_FromDouble(scalar);

		if (PyList_SetItem(temp, 0, py_scalar) == -1) /*steals the reference of py_scalar*/
		{
			return NULL;
		}

		PyTuple_SetItem(argList, 0, temp); /*steals the reference of temp*/

		array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);
	}
	else 
	{
		/*Passed*/
		array = (xparrayObject *)array_arg;

		/*Check if the multiplication can be done*/
		if ((self->rows != array->rows) || (self->cols != array->cols))
		{
			PyErr_SetString(PyExc_TypeError, "Both arguments must have the same size.");
			return NULL;
		}
	}

	/*Passed, initialize a new array with zeros*/
	PyObject *argList = PyTuple_New(self->rows);

	int i, j;
	for (i = 0; i < self->rows; i++)
	{
		PyObject *temp = PyList_New(self->cols);

		for (j = 0; j < self->cols; j++)
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
			ret_array->data[i][j] = self->data[i][j] + array->data[i*(1 - scalar_flag)][j*(1 - scalar_flag)];
		}
	}

	return ret_array;
}

xparrayObject *py_minus(xparrayObject *self, PyObject *args)
{
	PyObject *array_arg = NULL;
	double scalar;

	if (!PyArg_ParseTuple(args, "O", &array_arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(array_arg);

	if (PyXParray_Check(array_arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "The argument's type must be xnlpy.array or numeric.");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		array_arg = PyNumber_Float(array_arg);
		scalar = PyFloat_AsDouble(array_arg);

		/*Now create the array ([scalar])*/
		PyObject *argList = PyTuple_New(1);
		PyObject *temp = PyList_New(1);
		PyObject *py_scalar = PyFloat_FromDouble(scalar);

		if (PyList_SetItem(temp, 0, py_scalar) == -1) /*steals the reference of py_scalar*/
		{
			return NULL;
		}

		PyTuple_SetItem(argList, 0, temp); /*steals the reference of temp*/

		array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);
	}
	else 
	{
		/*Passed*/
		array = (xparrayObject *)array_arg;

		/*Check if the multiplication can be done*/
		if ((self->rows != array->rows) || (self->cols != array->cols))
		{
			PyErr_SetString(PyExc_TypeError, "Both arguments must have the same size.");
			return NULL;
		}
	}

	/*Passed, initialize a new array with zeros*/
	PyObject *argList = PyTuple_New(self->rows);

	int i, j;
	for (i = 0; i < self->rows; i++)
	{
		PyObject *temp = PyList_New(self->cols);

		for (j = 0; j < self->cols; j++)
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
			ret_array->data[i][j] = self->data[i][j] - array->data[i*(1 - scalar_flag)][j*(1 - scalar_flag)];
		}
	}

	return ret_array;
}

static PyMethodDef xparray_methods[] = 
{
	{"read", (PyCFunction) xparray_read, METH_VARARGS, "Return the element at (row, column)"},
	{"write", (PyCFunction) xparray_write, METH_VARARGS, "Set the element at (row, column)"},
	{"print", (PyCFunction) xparray_print, METH_VARARGS, "Print the array"},
	{"dot", (PyCFunction) py_dot, METH_VARARGS, "Array multiplication"},
	{"plus", (PyCFunction) py_plus, METH_VARARGS, "Array addition"},
	{"minus", (PyCFunction) py_minus, METH_VARARGS, "Array subtraction"},
	{NULL} /*Sentinel*/
};

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
};