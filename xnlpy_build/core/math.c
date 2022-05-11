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

#include <math.h>
#include <stdio.h>

/**
 * 
 * acos function
 * 
 */
PyObject *PyXP_Acos(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", acos(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = acos(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * acosg function
 * 
 */
PyObject *PyXP_Acosh(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", acosh(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = acosh(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * asin function
 * 
 */
PyObject *PyXP_Asin(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", asin(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = asin(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * asinh function
 * 
 */
PyObject *PyXP_Asinh(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", asinh(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = asinh(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * atan function
 * 
 */
PyObject *PyXP_Atan(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", atan(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = atan(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * atanh function
 * 
 */
PyObject *PyXP_Atanh(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", atanh(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = atanh(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * cos function
 * 
 */
PyObject *PyXP_Cos(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", cos(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = cos(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * cosh function
 * 
 */
PyObject *PyXP_Cosh(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", cosh(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = cosh(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * sin function
 * 
 */
PyObject *PyXP_Sin(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", sin(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = sin(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * sinh function
 * 
 */
PyObject *PyXP_Sinh(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", sinh(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = sinh(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * tan function
 * 
 */
PyObject *PyXP_Tan(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", tan(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = tan(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * tanh function
 * 
 */
PyObject *PyXP_Tanh(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", tanh(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = tanh(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * exp function
 * 
 */
PyObject *PyXP_Exp(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", exp(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = exp(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * log function
 * 
 */
PyObject *PyXP_Log(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", log(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = log(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * log10 function
 * 
 */
PyObject *PyXP_Log10(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", log10(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = log10(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * sqrt function
 * 
 */
PyObject *PyXP_Sqrt(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", sqrt(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = sqrt(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * ceil function
 * 
 */
PyObject *PyXP_Ceil(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", ceil(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = ceil(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * fabs function
 * 
 */
PyObject *PyXP_Fabs(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", fabs(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = fabs(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

/**
 * 
 * floor function
 * 
 */
PyObject *PyXP_Floor(PyObject *self, PyObject *args)
{
	double scalar;
	PyObject *arg = NULL;

	if (!PyArg_ParseTuple(args,"O", &arg))
	{
		return NULL;
	}

	int scalar_flag = PyNumber_Check(arg);

	if (PyXParray_Check(arg) && (scalar_flag != 1))
	{
		PyErr_SetString(PyExc_TypeError, "You must provide a xnl array or a number as argument.\n");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array = NULL;

	if (scalar_flag == 1)
	{
		scalar = PyFloat_AsDouble(PyNumber_Float(arg));
		
		return Py_BuildValue("d", floor(scalar));
	}

	/*xnl array as argument*/

	array = (xparrayObject *) arg;

	xparrayObject *ret_array = NULL;

	ret_array = PyXParray_New(array->rows, array->cols);

	Py_ssize_t i, j;

	for (i = 0; i < array->rows; i++)
	{
		for (j = 0; j < array->cols; j++)
		{
			ret_array->data[i][j] = floor(array->data[i][j]);
		}
	}

	return (PyObject *) ret_array;
}

