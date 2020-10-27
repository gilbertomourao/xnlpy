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
static double xnl_diff(double (*f)(double), double x, double h, int iterations)
{
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
		Rdiff[i][0] = (f(x + h) - f(x - h)) / (2 * h); /*Central Divided Differentiation*/
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

	return retval;
}

/*Top level function*/
static double diff(double (*f)(double), double x, double h, int iterations)
{
	/*Check if the input iterations is correct*/
	if (iterations <= 0)
	{
		printf("The number of iterations must be positive.\n");
		exit(EXIT_FAILURE);
	}

	return xnl_diff(f, x, h, iterations);
}

/**********************************************
* Python interface
***********************************************/

static PyObject *diff_cb_callable;

/*Implementation of double (*func)(double)*/

static double diff_callback(double callback)
{
	PyObject *retval;
	double result; /*returned value*/

	PyObject *value = Py_BuildValue("d", callback);
	PyObject *arglist = PyTuple_New(1);

	PyTuple_SetItem(arglist, 0, value); /*steals the reference of value*/

	/*call the python function/object saved below*/
	retval = PyObject_CallObject(diff_cb_callable, arglist);

	/*convert the returned object to a double, if possible*/
	if (retval && PyFloat_Check(retval))
	{
		result = PyFloat_AsDouble(retval);
	}
	else
	{
		PyErr_SetString(PyExc_TypeError, "Callable object: bad behavior");
		return -1; /* Error. You must propagate this.*/
	}

	Py_XDECREF(retval);
	Py_DECREF(arglist);

	return result;
}

/* function to be called externally */
PyObject *PyXP_Diff(PyObject *self, PyObject *args, PyObject *kwargs)
{
	/*arguments*/
	double x;
	double h = 0.01;
	int iterations = 5;
	/*values to be returned*/
	double result;

	static char *kwlist[] = {"f", "x", "step", "iterations", NULL};

	/* parse the input tuple with keywords */

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Od|$dI", kwlist, &diff_cb_callable, &x, &h, &iterations))
	{
		return NULL;
	}

	/* ensure the first parameter is a callable */
	if (!PyCallable_Check(diff_cb_callable))
	{
		PyErr_SetString(PyExc_TypeError, "diff: a callable is required\n");
		return NULL;
	}

	/* Call the diff function */
	result = diff(diff_callback,x,h,iterations);

	return Py_BuildValue("d", result);
}
