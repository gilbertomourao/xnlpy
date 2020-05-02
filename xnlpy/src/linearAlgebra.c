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

xparrayObject *py_zeros(PyObject *self, PyObject *args)
{
	int rows, cols;
	if (!PyArg_ParseTuple(args, "ii", &rows, &cols))
	{
		return NULL;
	}

	if (rows <= 0 || cols <= 0)
	{
		PyErr_SetString(PyExc_TypeError, "the arguments must be positive integers.");
		return NULL;
	}

	PyObject *argList = PyTuple_New(rows);

	int i;
	Py_ssize_t j, ListSize = cols;
	for (i = 0; i < rows;i++)
	{
		PyObject *temp = PyList_New(ListSize);

		for (j = 0; j < ListSize; j++)
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
	xparrayObject *array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (array == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	/*Success on the creation, then return*/

	return array;
}

xparrayObject *py_ones(PyObject *self, PyObject *args)
{
	int rows, cols;
	if (!PyArg_ParseTuple(args, "ii", &rows, &cols))
	{
		return NULL;
	}

	if (rows <= 0 || cols <= 0)
	{
		PyErr_SetString(PyExc_TypeError, "the arguments must be positive integers.");
		return NULL;
	}

	PyObject *argList = PyTuple_New(rows);

	int i;
	Py_ssize_t j, ListSize = cols;
	for (i = 0; i < rows;i++)
	{
		PyObject *temp = PyList_New(ListSize);

		for (j = 0; j < ListSize; j++)
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
	xparrayObject *array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (array == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	/*Success on the creation, then populate with 1*/

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			array->data[i][j] = 1;
		}
	}

	return array;
}

xparrayObject *py_eye(PyObject *self, PyObject *args)
{
	int size;
	if (!PyArg_ParseTuple(args, "i", &size))
	{
		return NULL;
	}

	if (size <= 0)
	{
		PyErr_SetString(PyExc_TypeError, "the argument must be a positive integer.");
		return NULL;
	}

	PyObject *argList = PyTuple_New(size);

	int i;
	Py_ssize_t j, ListSize = size;
	for (i = 0; i < size;i++)
	{
		PyObject *temp = PyList_New(ListSize);

		for (j = 0; j < ListSize; j++)
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
	xparrayObject *array = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (array == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	/*Success on the creation, then populate with 1 where i == j*/

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			array->data[i][j] = (i == j);
		}
	}

	return array;
}

xparrayObject *py_mult(PyObject *self, PyObject *args)
{
	PyObject *arg_A = NULL;
	PyObject *arg_B = NULL;

	if (!PyArg_ParseTuple(args, "OO", &arg_A, &arg_B))
	{
		return NULL;
	}

	if (PyXParray_Check(arg_A) || PyXParray_Check(arg_B))
	{
		PyErr_SetString(PyExc_TypeError, "The arguments' type must be xnlpy.array.");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array_A = (xparrayObject *)arg_A;
	xparrayObject *array_B = (xparrayObject *)arg_B;

	/*Check if the multiplication can be done*/
	if (array_A->cols != array_B->rows)
	{
		PyErr_SetString(PyExc_TypeError, "The number of columns of A must be equal to the number of rows of B.");
		return NULL;
	}

	/*Passed, initialize a new array with zeros*/
	PyObject *argList = PyTuple_New(array_A->rows);

	int i, j;
	for (i = 0; i < array_A->rows; i++)
	{
		PyObject *temp = PyList_New(array_B->cols);

		for (j = 0; j < array_B->cols; j++)
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
	xparrayObject *array_C = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (array_C == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	/*Success on the creation, then do the multiplication*/
	int k;

	for (i = 0; i < array_A->rows; i++)
	{
		for (j = 0; j < array_B->cols; j++)
		{
			for (k = 0; k < array_B->rows; k++)
			{
				array_C->data[i][j] += array_A->data[i][k] * array_B->data[k][j];
			}
		}
	}

	return array_C;
}

xparrayObject *py_transpose(PyObject *self, PyObject *args)
{
	PyObject *arg_A = NULL;

	if (!PyArg_ParseTuple(args, "O", &arg_A))
	{
		return NULL;
	}

	if (PyXParray_Check(arg_A))
	{
		PyErr_SetString(PyExc_TypeError, "The argument type must be xnlpy.array.");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array_A = (xparrayObject *)arg_A;

	/*Initialize a new array with zeros*/
	PyObject *argList = PyTuple_New(array_A->cols);

	int i, j;
	for (i = 0; i < array_A->cols; i++)
	{
		PyObject *temp = PyList_New(array_A->rows);

		for (j = 0; j < array_A->rows; j++)
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
	xparrayObject *array_Atr = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (array_Atr == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	/*Success on the creation, then return the transpose array*/
	for (i = 0; i < array_Atr->rows; i++)
	{
		for (j = 0; j < array_Atr->cols; j++)
		{
			array_Atr->data[i][j] = array_A->data[j][i];
		}
	}

	return array_Atr;
}

xparrayObject *py_add(PyObject *self, PyObject *args)
{
	PyObject *arg_A = NULL;
	PyObject *arg_B = NULL;

	if (!PyArg_ParseTuple(args, "OO", &arg_A, &arg_B))
	{
		return NULL;
	}

	if (PyXParray_Check(arg_A) || PyXParray_Check(arg_B))
	{
		PyErr_SetString(PyExc_TypeError, "The arguments' type must be xnlpy.array.");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array_A = (xparrayObject *)arg_A;
	xparrayObject *array_B = (xparrayObject *)arg_B;

	/*Check if the multiplication can be done*/
	if ((array_A->rows != array_B->rows) || (array_A->cols != array_B->cols))
	{
		PyErr_SetString(PyExc_TypeError, "Both arguments must have the same size.");
		return NULL;
	}

	/*Passed, initialize a new array with zeros*/
	PyObject *argList = PyTuple_New(array_A->rows);

	int i, j;
	for (i = 0; i < array_A->rows; i++)
	{
		PyObject *temp = PyList_New(array_B->cols);

		for (j = 0; j < array_B->cols; j++)
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
	xparrayObject *array_C = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (array_C == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	/*Success on the creation, then do the addition*/
	for (i = 0; i < array_C->rows; i++)
	{
		for (j = 0; j < array_C->cols; j++)
		{
			array_C->data[i][j] = array_A->data[i][j] + array_B->data[i][j];
		}
	}

	return array_C;
}

xparrayObject *py_sub(PyObject *self, PyObject *args)
{
	PyObject *arg_A = NULL;
	PyObject *arg_B = NULL;

	if (!PyArg_ParseTuple(args, "OO", &arg_A, &arg_B))
	{
		return NULL;
	}

	if (PyXParray_Check(arg_A) || PyXParray_Check(arg_B))
	{
		PyErr_SetString(PyExc_TypeError, "The arguments' type must be xnlpy.array.");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array_A = (xparrayObject *)arg_A;
	xparrayObject *array_B = (xparrayObject *)arg_B;

	/*Check if the multiplication can be done*/
	if ((array_A->rows != array_B->rows) || (array_A->cols != array_B->cols))
	{
		PyErr_SetString(PyExc_TypeError, "Both arguments must have the same size.");
		return NULL;
	}

	/*Passed, initialize a new array with zeros*/
	PyObject *argList = PyTuple_New(array_A->rows);

	int i, j;
	for (i = 0; i < array_A->rows; i++)
	{
		PyObject *temp = PyList_New(array_B->cols);

		for (j = 0; j < array_B->cols; j++)
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
	xparrayObject *array_C = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

	Py_DECREF(argList);

	if (array_C == NULL)
	{
		/*Will print the messages from xparray*/
		return NULL;
	}

	/*Success on the creation, then do the addition*/
	for (i = 0; i < array_C->rows; i++)
	{
		for (j = 0; j < array_C->cols; j++)
		{
			array_C->data[i][j] = array_A->data[i][j] - array_B->data[i][j];
		}
	}

	return array_C;
}