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
#include <math.h>
/**
 * Numerical differentiation by Richardson Extrapolation
 *
 * Evaluate the first derivative of the function at x
 */
static double xnl_diff(double (*f)(xparrayObject *, PyObject *), PyObject *f_args, double x, double h, int iterations)
{
	xparrayObject *arg_array = NULL;

	arg_array = PyXParray_New(1, 1);

	/* Set the first element of the tuple f_args to be the array argument (x) */
	/* Insert a referente to arg_array at the position 0 */
	if (f_args)
		PyTuple_SetItem(f_args, 0, (PyObject *) arg_array); /* steals the reference of arg_array */

	int i, j;
	double retval;
	double **Rdiff = malloc(iterations * sizeof(double));

	/*Check if malloc worked*/
	if (Rdiff == NULL)
	{
		printf("ERROR: In xnl_diff. Memory allocation failed! Check variable Rdiff.\n");
		exit(EXIT_FAILURE);
	}

	/*Initialize Richardson matrix*/
	for (i = 0; i < iterations; i++)
	{
		Rdiff[i] = malloc(iterations * sizeof(double));
		/*Check if malloc worked*/
		if (Rdiff[i] == NULL)
		{
			printf("ERROR: In xnl_diff. Memory allocation failed! Check variable Rdiff[%d].\n", i);
			exit(EXIT_FAILURE);
		}
	}

	/*Richardson Extrapolation*/
	for (i = 0; i < iterations; i++)
	{
		/* Set function + side */
		arg_array->data[0][0] = x+h;
		double f_plus = f(arg_array, f_args);
		/* Set function - side */
		arg_array->data[0][0] = x-h;
		double f_minus = f(arg_array, f_args);

		Rdiff[i][0] = (f_plus - f_minus) / (2 * h); /*Central Divided Differentiation*/
		for (j = 1; j <= i; j++)
		{
			Rdiff[i][j] = (pow(4,j)*Rdiff[i][j-1] - Rdiff[i-1][j-1]) / (pow(4,j) - 1);
		}
		h /= 2;
	}

	retval = Rdiff[iterations - 1][iterations - 1];

	/*free memory*/
	for (i = 0; i < iterations; i++)
	{
		free(Rdiff[i]);
	}

	free(Rdiff);

	if (!f_args)
		Py_XDECREF(arg_array);

	return retval;
}

/*Top level function*/
static double diff(double (*f)(xparrayObject *, PyObject *), PyObject *f_args, double x, double h, int iterations)
{
	/* Handles error of args parameter */	
	if (PyErr_Occurred())
	{
		return 0; /* just returning anything */
	}

	/*Checks if the input iterations is correct*/
	if (iterations <= 0)
	{
		/* Set the error flag */
		PyErr_SetString(PyExc_TypeError, "The number of iterations must be positive.");
		return 0; /* just returning anything */
	}

	return xnl_diff(f, f_args, x, h, iterations);
}

/**********************************************
* Python interface
***********************************************/

static PyObject *diff_cb_callable = NULL;

/*Implementation of double (*func)(xparrayObject, PyObject *)*/

static double diff_callback(xparrayObject *arg_array, PyObject *f_args)
{
	PyObject *retval;
	double result; /*returned value*/

	/*calls the python function/object saved below*/
	if (!f_args)
	{
		retval = PyObject_CallFunctionObjArgs(diff_cb_callable, (PyObject *) arg_array, NULL);
	}
	else
	{
		retval = PyObject_CallObject(diff_cb_callable, f_args);
	}

	/*converts the returned object to a xparray, if possible*/
	if (retval && !PyXParray_Check(retval))
	{
		xparrayObject *ret_array = (xparrayObject *) retval;

		/* Not vectorized yet */
		result = ret_array->data[0][0];
	}
	else
	{
		/* sets the error flag */
		PyErr_SetString(PyExc_TypeError, "Callable object: bad behavior");
		return 0; /* Error. You must propagate this.*/
	}

	Py_XDECREF(retval);

	return result;
}

/**
 * Check if the args provided by the user are ok.
 */
static PyObject *check_args(PyObject *user_data)
{
	/**
	 * Verifies if user_data must be used. This is the case 
	 * when the function have more than one argument.
	 *
	 * To count the number of the callable arguments, refer to 
	 * https://stackoverflow.com/questions/1117164/how-to-find-the-number-of-parameters-to-a-python-function-from-c
	 */
	unsigned count = 0;
	PyObject *fc = PyObject_GetAttrString(diff_cb_callable, "__code__");
	if (fc)
	{
		PyObject *ac = PyObject_GetAttrString(fc, "co_argcount");
		if (ac)
		{
			count = PyLong_AsLong(ac);
			Py_DECREF(ac);
		}
		Py_DECREF(fc);
	}

	if (count == 1 && user_data)
	{
		PyErr_SetString(PyExc_TypeError, 
						"diff: args not available for this function. "
						"It must have more than one argument.");
		return NULL;
	}

	if (count > 1 && !user_data)
	{
		PyErr_SetString(PyExc_TypeError, 
						"diff: invalid number of arguments. The function "
						"appears to have more than one argument, so you must "
						"use \"args\" to pass the other arguments to it.");
		return NULL;
	}

	if (!user_data)
		return NULL;

	/* ensure the parameter args is a tuple */
	if (!PyTuple_Check(user_data))
	{
		PyErr_SetString(PyExc_TypeError, "diff: args must be a tuple");
		return NULL;
	}

	/* ensure args is a tuple of numbers */
	Py_ssize_t tuplen = PyTuple_Size(user_data);

	if (tuplen != (count - 1))
	{
		PyErr_SetString(PyExc_RuntimeError, "diff: args' size must be equal to "
											"the number of extra arguments of the "
											"function.");
		return NULL;
	}

	Py_ssize_t i;

	for (i = 0; i < tuplen; i++)
	{
		if (!PyNumber_Check(PyTuple_GetItem(user_data, i)))
		{
			PyErr_SetString(PyExc_TypeError, "diff: args must be a tuple of numbers");
			return NULL;
		}
	}

	/* user_data is ok. Now creates a tuple with the callable args */
	/**
	 * The first argument of f is a double. So it's size is 1 + tuplen.
	 */
	PyObject *f_args;

	f_args = PyTuple_New(1 + tuplen);

	for (i = 1; i <= tuplen; i++)
	{
		/* PyTuple_GetItem borrows the reference of user_data[i-1] */
		/* PyTuple_SetItem steals the reference of data */
		PyObject *data = PyNumber_Float(PyTuple_GetItem(user_data, i-1));
		PyTuple_SetItem(f_args, i, data);
	}

	return f_args;
}

/* function to be called externally */
PyObject *PyXP_Diff(PyObject *self, PyObject *args, PyObject *kwargs)
{
	/*arguments*/
	double x;
	PyObject *user_data = NULL;
	double h = 0.01;
	int iterations = 5;
	/*values to be returned*/
	double result;

	static char *kwlist[] = {"f", "x", "args", "step", "iterations", NULL};

	/* parse the input tuple with keywords */

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Od|$OdI", kwlist, &diff_cb_callable, &x, &user_data, &h, &iterations))
	{
		PyErr_SetString(PyExc_TypeError, "diff: failed to parse arguments");
		return NULL;
	}

	/* ensure the first parameter is a callable */
	if (!PyCallable_Check(diff_cb_callable))
	{
		PyErr_SetString(PyExc_TypeError, "diff: a callable is required");
		return NULL;
	}

	PyObject *f_args = check_args(user_data);

	/* Calls the diff function */
	result = diff(diff_callback,f_args,x,h,iterations);

	Py_XDECREF(f_args);

	/* Error occurred when working with diff.
	 * The only possibility is error due to a bad 
	 * behavior of the callable object. Just returns 
	 * NULL to indicate to python that an error occurred 
	 * without closing the interpreter.
	 */ 
	if (PyErr_Occurred())
	{
		return NULL;
	}

	/* Everything ok */
	return Py_BuildValue("d", result);
}
