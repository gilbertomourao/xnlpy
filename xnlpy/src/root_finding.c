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

/*Error handler*/
#define NO_DOUBLE_ARG 0.0f / 0.0f

/**
 * Brent's Method
 *
 * 
 */
static double brent_method(double (*f)(double), double x0, double x1, double tolerance, int max_iter)
{
	if (x0 != x0 || x1 != x1)
	{
		printf("ERROR: In brent_method. No argument passed to x0 or x1. You must pass arguments for both x0 and x1.\n");
		return NO_DOUBLE_ARG;
	}

	int counter;

	double a = x0;
	double b = x1;
	double fa = f(a);
	double fb = f(b);
	double fs = 0; /* initialize */

	if ( !(fa * fb < 0) )
	{
		printf("ERROR: In brent_method. Signs of f(lower_bound) and f(upper_bound) must be opposites.\n");
		return NO_DOUBLE_ARG;
	}

	if (fabs(fa) < fabs(fb))
	{
		double aux = a;

		a = b;
		b = aux;

		aux = fa;
		fa = fb;
		fb = aux;
	}

	double c = a;
	double fc = fa;
	int mflag = 1;
	double s = 0;
	double d = 0;

	for (counter = 1; counter < max_iter; counter++)
	{
		/*check tolerance*/
		if (fabs(b-a) < tolerance)
		{
			return s;
		}

		if (fa != fc && fb != fc)
		{
			/*inverse quadratic interpolation*/
			s = ( a * fb * fc / ((fa - fb) * (fa - fc)) )
				+ ( b * fa * fc / ((fb - fa) * (fb - fc)) )
				+ ( c * fa * fb / ((fc - fa) * (fc - fb)) );
		}
		else
		{
			/*secant method*/
			s = b - fb * (b - a) / (fb - fa);
		}

		/* checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection */
		if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
				( mflag && (fabs(s-b) >= (fabs(b-c) * 0.5)) ) ||
				( !mflag && (fabs(s-b) >= (fabs(c-d) * 0.5)) ) ||
				( mflag && (fabs(b-c) < tolerance) ) ||
				( !mflag && (fabs(c-d) < tolerance))	)
		{
			/* bisection method */
			s = (a+b)*0.5;
 
			mflag = 1;
		}
		else
		{
			mflag = 0;
		}

		fs = f(s);
		d = c;
		c = b;
		fc = fb;

		if (fa * fs < 0)
		{
			b = s;
			fb = fs;
		}
		else 
		{
			a = s;
			fa = fs;
		}

		if (fabs(fa) < fabs(fb))
		{
			double aux = a;

			a = b;
			b = aux;

			aux = fa;
			fa = fb;
			fb = aux;
		}
	}

	printf("WARNING: In brent_method. Reached the maximum number of iterations.\n");

	return s;
}

/**
 * Newton-Raphson Method
 *
 * 
 */
static double newton_method(double (*f)(double), double (*df)(double), double x0, double tolerance, int max_iter)
{
	if (x0 != x0)
	{
		printf("ERROR: In newton_method. No argument passed to x0. You must pass an argument for x0 to use this function.\n");
		return NO_DOUBLE_ARG;
	}

	double h, x1;
	int counter = 0;

	do
	{
		h = f(x0) / df(x0);
		x1 = x0 - h;
		x0 = x1;
		++counter;
	} while (fabs(h) >= tolerance && counter < max_iter);

	if (counter == max_iter)
	{
		printf("WARNING: In newton_method. Reached the maximum number of iterations.\n");
	}

	return x1;
}

/**
 * Kernel of fsolve function
 *
 * 
 */
static double kernel_fsolve(double (*f)(double), double (*df)(double), double x0, double x1, double tolerance, int max_iter)
{
	if (f == NULL)
	{
		printf("ERROR: In kernel_fsolve. You must pass a function to the f argument.\n");
		return NO_DOUBLE_ARG;
	}
	return (df == NULL) ? brent_method(f, x0, x1, tolerance, max_iter) : newton_method(f, df, x0, tolerance, max_iter);
}

/**********************************************
* Python interface
***********************************************/

/*Implementation of double (*func)(double)*/
static PyObject *f_cb_callable = NULL;

static double f_callback(double callback)
{
	PyObject *retval;
	double result; /*returned value*/

	PyObject *value = Py_BuildValue("d", callback);
	PyObject *arglist = PyTuple_New(1);

	PyTuple_SetItem(arglist, 0, value); /*steals the reference of value*/

	/*call the python function/object saved below*/
	retval = PyObject_CallObject(f_cb_callable, arglist);

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

/*Implementation of double (*dfunc)(double)*/
static PyObject *df_cb_callable = NULL;

static double df_callback(double callback)
{
	PyObject *retval;
	double result; /*returned value*/

	PyObject *value = Py_BuildValue("d", callback);
	PyObject *arglist = PyTuple_New(1);

	PyTuple_SetItem(arglist, 0, value); /*steals the reference of value*/

	/*call the python function/object saved below*/
	retval = PyObject_CallObject(df_cb_callable, arglist);

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

PyObject *PyXP_Fsolve(PyObject *self, PyObject *args, PyObject *kwargs)
{
	/*arguments*/
	double x0 = NO_DOUBLE_ARG;
	double x1 = NO_DOUBLE_ARG;
	PyObject *x_list;
	double tolerance = 1e-6;
	int max_iter = 100;
	/*values to be returned*/
	double root;

	static char *kwlist[] = {"f", "df", "x", "tolerance", NULL};

	/* parse the input tuple with keywords */

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|$OOdI", kwlist, &f_cb_callable, &df_cb_callable, &x_list, &tolerance, &max_iter))
	{
		return NULL;
	}

	/* ensure the f parameter is a callable */
	if (!PyCallable_Check(f_cb_callable))
	{
		PyErr_SetString(PyExc_TypeError, "fsolve: a callable is required\n");
		return NULL;
	}

	/*checking if the x argument is a list*/
	if (!PyList_CheckExact(x_list))
	{
		PyErr_SetString(PyExc_TypeError, "The x argument must be a list.\n");
		return NULL;
	}

	Py_ssize_t list_size = PyList_Size(x_list);

	if (df_cb_callable == NULL)
	{
		if (list_size != 2)
		{
			PyErr_SetString(PyExc_TypeError, "fsolve requires two arguments on the x list when the jacobian is not provided.\n");
			return NULL;
		}

		PyObject *py_x0 = PyList_GetItem(x_list, 0); /*borrowed reference to returned value*/
		PyObject *py_x1 = PyList_GetItem(x_list, 1);
		/*Check if list[0] is numeric*/
		if (PyNumber_Check(py_x0) != 1 || PyNumber_Check(py_x1) != 1)
		{
			printf("ERROR: In fsolve. The x list argument must have only numeric elements.\n");
			return NULL;
		}

		/*Convert list[0] and list[1] to python floats*/
		py_x0 = PyNumber_Float(py_x0);
		py_x1 = PyNumber_Float(py_x1);
		/*now to C double*/
		x0 = PyFloat_AsDouble(py_x0);
		x1 = PyFloat_AsDouble(py_x1);

		root = kernel_fsolve(f_callback, NULL, x0, x1, tolerance, max_iter);
	}
	else
	{
		/* ensure the df parameter is a callable */
		if (!PyCallable_Check(df_cb_callable))
		{
			PyErr_SetString(PyExc_TypeError, "fsolve: a callable is required\n");
			return NULL;
		}	

		if (list_size != 1)
		{
			PyErr_SetString(PyExc_TypeError, "fsolve requires one argument on the x list when the jacobian is provided.\n");
			return NULL;
		}		

		PyObject *py_x0 = PyList_GetItem(x_list, 0);
		py_x0 = PyNumber_Float(py_x0);
		x0 = PyFloat_AsDouble(py_x0);

		/* Call the fsolve function */
		root = kernel_fsolve(f_callback, df_callback, x0, x1, tolerance, max_iter);
	}

	return Py_BuildValue("d", root);

}