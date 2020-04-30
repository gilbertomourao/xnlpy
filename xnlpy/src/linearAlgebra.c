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