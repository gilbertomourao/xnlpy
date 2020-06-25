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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define XNL_INFINITY HUGE_VAL
#define XNL_NAN (0.0f / 0.0f)

/**********************************************
* MATLAB array functions
***********************************************/

static double sum(double *vec, int vec_size)
{
	double vec_sum = 0;
	int i;

	for (i = 0;i < vec_size;i++)
	{
		vec_sum += vec[i];
	}

	return vec_sum;
}

static double inner_prod(double *vec1, double *vec2, int vec_size)
{
	double prod = 0;
	unsigned i;

	for (i = 0;i < vec_size;i++)
	{
		prod += vec1[i] * vec2[i];
	}

	return prod;
}

static void cumprod(double value, double *vec, unsigned begin, unsigned end)
{
	unsigned i;

	vec[begin] = value;

	for (i = begin + 1;i <= end;i++)
	{
		vec[i] = value * vec[i-1];
	}
}

static void reverse(double *vec, int vec_size)
{
	unsigned i;
	double aux;

	for (i = 0;i < vec_size / 2; i++)
	{
		aux = vec[i];
		vec[i] = vec[vec_size - 1 - i];
		vec[vec_size - 1 - i] = aux;
	}
}

/**********************************************
* Legendre points using GLR algorithm
***********************************************/

static void GLR_compute_start(int n, double *p, double *pp)
{
	double pm1 = 1.0;
	double pm2 = 0.0;
	double ppm1 = 0.0;
	double ppm2 = 0.0;
	double dk;
	int k;

	for (k = 0; k < n;k++)
	{
		dk = (double) k;
		*p = - dk * pm2 / (dk + 1);
		*pp = ( (2 * dk + 1) * pm1 - dk * ppm2 ) / (dk + 1);
		pm2 = pm1;
		pm1 = *p;
		ppm2 = ppm1;
		ppm1 = *pp;
	}
}

static void GLR_compute_first_root(double *p, double *roots, double *weights, int n)
{
	const double pi = 3.14159265358979323846;
	const double eps = 2.220446049250313e-16;

	/* Francesco Tricomi formula */
	int k = (n / 2) + (n % 2);
	double dk = (double) k;
	double dn = (double) n;
	double theta = pi * (4 * dk - 1.0) / (4 * dn + 2.0);
	double sine = sin(theta);
	double cosine = cos(theta);
	double x1 = (1.0 - (dn - 1) / (8*dn*dn*dn) - 1.0 / (384*dn*dn*dn*dn) * (39.0 - 28.0 / (sine * sine))) * cosine;

	/* Number of terms in Taylor expansion */
	int m = 30;
	double M = 1.0 / x1; /* scaling */

	/* For Newton Iteration */
	double step = 1.0;
	int l = 0;

	/* Recurrence relation for Legendre Polynomials */
	double *u = malloc((m+1) * sizeof(double));
	double *up = malloc((m+1) * sizeof(double));
	double *x1k = malloc((m+1) * sizeof(double));

	/*Check if malloc worked*/
	if (u == NULL)
	{
		printf("ERROR: In GLR_compute_first_root. Memory allocation failed! Check variable u.\n");
		exit(EXIT_FAILURE);
	}

	if (up == NULL)
	{
		printf("ERROR: In GLR_compute_first_root. Memory allocation failed! Check variable up.\n");
		exit(EXIT_FAILURE);
	}

	if (x1k == NULL)
	{
		printf("ERROR: In GLR_compute_first_root. Memory allocation failed! Check variable x1k.\n");
		exit(EXIT_FAILURE);
	}

	/*Initialize u and up with 0 and xk with 1*/

	for (k = 0;k < m + 1;k++)
	{
		u[k] = up[k] = 0;
		x1k[k] = 1;
	}

	u[0] = *p;

	for (k = 0;k < m - 1;k += 2)
	{
		dk = (double) k;
		u[k + 2] = (dk - dn * (dn + 1) / (dk + 1)) * u[k] / (M * M * (dk + 2));
		up[k + 1] = (dk + 2) * u[k + 2] * M;
	}

	/* Flip for more accuracy in inner product calculation */

	reverse(u, m+1);
	reverse(up, m+1);

	/* Newton iteration */

	while (eps < fabs(step) && l < 10)
	{
		l += 1;
		step = inner_prod(u,x1k,m+1) / inner_prod(up,x1k,m+1);
		x1 -= step;
		x1k[0] = 1.0;
		cumprod(M * x1, x1k, 1, m);
		reverse(x1k, m+1);
	}

	roots[n/2] = x1;
	weights[n/2] = inner_prod(up,x1k,m+1);

	/* Free memory */

	free(u);
	free(up);
	free(x1k);
}

static void GLR_compute_all(double *roots, double *weights, int n)
{
	const double pi = 3.14159265358979323846;
	const double eps = 2.220446049250313e-16;

	int N;
	int s;
	int k;

	/* Francesco tricomi formula */

	double theta;
	double sine;
	double cosine;
	double dn = (double) n;
	double dk;

	/* Number of terms in Taylor expansion */

	int m = 30;
	double *u = malloc((m + 1) * sizeof(double));
	double *up = malloc((m + 1) * sizeof(double));
	double *hh = malloc((m + 1) * sizeof(double));

	/*Check if malloc worked*/
	if (u == NULL)
	{
		printf("ERROR: In GLR_compute_all. Memory allocation failed! Check variable u.\n");
		exit(EXIT_FAILURE);
	}

	if (up == NULL)
	{
		printf("ERROR: In GLR_compute_all. Memory allocation failed! Check variable up.\n");
		exit(EXIT_FAILURE);
	}

	if (hh == NULL)
	{
		printf("ERROR: In GLR_compute_all. Memory allocation failed! Check variable hh.\n");
		exit(EXIT_FAILURE);
	}

	double xp;
	double h;
	double M;

	/* For the final loop */

	int j;
	int it_h;

	/* For Newton's algorithm*/

	double step;
	int l;

	/* Check if n is odd or even */

	if (n % 2)
	{
		N = (n - 1) / 2;
		s = 1;
	}
	else 
	{
		N = n / 2;
		s = 0;
	}

	/* Tricomi's algorithm */

	for ( k = (n-2+s)/2 ; k ; k-- )
	{
		dk = (double) k;
		theta = pi * (4 * dk - 1.0) / (4 * dn + 2.0);
		sine = sin(theta);
		cosine = cos(theta);
		roots[n - k] = (1.0 - (dn - 1.0) / (8*dn*dn*dn) - 1.0 / (384*dn*dn*dn*dn) * (39.0 - 28.0 / (sine * sine))) * cosine;
	}

	/*Initialize u and up with 0*/
	
	for (k = 0;k < m + 1;k++)
	{
		u[k] = up[k] = 0;
	}

	/* Loop algorithm */

	for (j = N;j < n - 1;j++)
	{
		/* revious root */
		xp = roots[j];

		/* Initial approximation via Tricomi's formula */
		h = roots[j+1] - xp;

		/* scaling */
		M = 1.0 / h;

		/* recurrence relation for Legendre Polynomials */
		u[0] = 0.0;
		u[1] = weights[j] / M;
		up[0] = u[1];
		up[m] = 0;

		for (k = 0;k < m - 1;k++)
		{
			dk = (double) k;
			u[k+2] = (2 * xp * (dk + 1)/M * u[k+1] + (dk - dn*(dn+1)/(dk+1)) * u[k] / M / M)
					/ ((1.0 - xp * xp) * (dk + 2));
			up[k+1] = (dk+2) * u[k+2] * M;		
		}

		/* Flip for more accuracy in inner product calculation */
		reverse(u,m+1);
		reverse(up,m+1);

		/* Initialize h with 1 */

		for (it_h = 0;it_h < m + 1;it_h++)
		{
			hh[it_h] = 1;
		}

		hh[m] = M;

		/* Newton iteration */

		step = 1.0;
		l = 0;

		while (eps < fabs(step) && l < 10)
		{
			l += 1;
			step = inner_prod(u,hh,m+1) / inner_prod(up,hh,m+1);
			h -= step;
			hh[0] = M;
			cumprod(M * h, hh, 1, m);
			reverse(hh,m+1);
		}

		/* Update */
		roots[j+1] = xp + h;
		weights[j+1] = inner_prod(up,hh,m+1);
	}

	/* Nodes are symmetric */

	for (k = 0;k < N + s;k++)
	{
		roots[k] = - roots[n-k-1];
		weights[k] = weights[n-k-1];
	}

	/* Free memory */

	free(u);
	free(up);
	free(hh);
}

static void GLR_compute_points(double *roots, double *weights, int n)
{
	double p = 0;
	double pp = 0;
	int i; /*for loops at the end*/
	double w_sum;

	GLR_compute_start(n, &p, &pp);

	if (n % 2) /* odd */
	{
		roots[(n + 1)/2 - 1] = 0.0;
		weights[(n + 1)/2 - 1] = pp;
	}
	else 
	{
		GLR_compute_first_root(&p, roots, weights, n);
	}

	GLR_compute_all(roots, weights, n);

	for (i = 0;i < n;i++)
	{
		weights[i] = 2.0 / (1.0 - roots[i]) / (1.0 + roots[i]) / weights[i] / weights[i];
	}

	w_sum = sum(weights, n);

	for (i = 0;i < n;i++)
	{
		weights[i] *= 2.0 / w_sum;
	}
}

/*Legendre Polynomial and GLR algorithm functions*/

static void create(double *roots, double *weights, int Order)
{
	int i;

	/*Initialize the arrays with 0*/
	/*The GLR algorithm presented here depends on this initialization*/
	for (i = 0;i < Order;i++)
	{
		roots[i] = weights[i] = 0;
	}

	/*Apply the GLR algorithm*/
	GLR_compute_points(roots, weights, Order);
}

/************************************************
* Function that calculates the approximated value
* of the integral using Gauss-Legendre quadrature
*************************************************/

static double GaussLegendre(double (*f)(double), double a, double b, int points, double *roots, double *weights)
{
	double c1 = (b - a)/2, c2 = (b + a)/2;
    double total = 0;
    int i;

    for (i = 0; i < points; ++i) 
    {
        total += weights[i] * f(c1 * roots[i] + c2);
    }
    
    return c1 * total;
}

/**********************************************
* Function that applies the adaptive quadrature
* method
***********************************************/

static int src_isinf(double x)
{
    if (x == XNL_INFINITY || x == -XNL_INFINITY)
        return 1;
    else
        return 0;   
}

static int src_isnan(double x)
{
    return x != x;
}

static double Adaptive_Quad(double (*f)(double), /* function to integrate */
                   	 		double a, /* left interval boundary */
                   	 		double b, /* right interval boundary */
							int points, /* number of points for interpolation */
                   	 		double total, /* total interval integral value */
                   	 		double tolerance, /* error tolerance */
                   	 		int remaining_it, /* maximum number of iterations per division */
                   	 		double* error, /* stores the error value */
							double *roots, /* Legendre Polynomial roots */
							double *weights) /* Quadrature weights */
{
	double retvalL = 0, retvalR = 0;
	double left, right, int_error;

	/* Checks if there is remaining iterations */
	if (!remaining_it)
	{
	    *error = XNL_NAN;

	    return XNL_NAN;
	}

	left = GaussLegendre(f, a, (a+b)/2, points, roots, weights); /* integrate over the left interval */
	right = GaussLegendre(f, (a+b)/2, b, points, roots, weights); /* integrate over the right interval */

	int_error = fabs(total - (left + right)); /* interval error */

	if (int_error < tolerance) 
	{
		/* Satisfies the tolerance condition */
	    *error += int_error;

	    return left + right;
	}
	else 
	{
		/* Computes the Left node */
	    retvalL = Adaptive_Quad(f, a, (a + b)/2, points, left, tolerance, remaining_it - 1, error, roots, weights);

	    if (src_isnan(retvalL))
	    {
	        /* Checks if the left node isn't ok */
	        *error = XNL_NAN;
	    
	    	return retvalL;
	    } 
	    else 
	    {
	    	/* The left node is ok. Then the right node can be computed */
	        retvalR = Adaptive_Quad(f, (a + b)/2, b, points, right, tolerance, remaining_it - 1, error, roots, weights);

	    	return retvalL + retvalR; /* returns the left + right nodes */
	    }
	}
}

/**********************************************
* Adaptive gaussian quadrature
***********************************************/

static double (*func_pointer)(double);

static double g(double x)
{
	return func_pointer(tan(x))/(cos(x) * cos(x));
}

/*
 * Function: xnl_integral
 * Evaluates the integral of a given function f
 * on a given interval (a,b) with default
 * tolerance.
 * Parameters:
 *   f - integrand
 *   a - left interval boundary (can be +-INFINITY)
 *   b - right interval boundary (can be +-INFINITY)
 *   tolerance - error tolerance with default value of 0
 *   remaining_it - maximum number of divisions per subinterval
 *   *error - a pointer to a variable that will store the estimated
 *            error value. In general, this error is much bigger
 *            than the real one.
 * This function can handle infinite intervals with
 * the trigonometric substitution x = tan(u). If the
 * integrand possess a singularity on the interval,
 * the quadrature may not work well.
 */

static double xnl_integral(double (*f)(double), /* function to integrate */
               	 		   double a, /* left interval boundary */
               			   double b, /* right interval boundary */
						   int points, /* number of points for interpolation */
               			   double tolerance, /* error tolerance */
               			   int remaining_it, /* maximum number of iterations per division */
               			   double* error) /* stores the error value */
{
	const double pi = 3.14159265358979323846;
    double retval; /* Control variable */

    double est_error = 0; /* Estimated error */

    /* Legendre Polynomial roots and weights */

    double *roots = malloc(points * sizeof(double));
    double *weights = malloc(points * sizeof(double));

	/*Check if malloc worked*/
	if (roots == NULL)
	{
		printf("ERROR: In xnl_integral. Memory allocation failed! Check variable roots.\n");
		exit(EXIT_FAILURE);
	}

	if (weights == NULL)
	{
		printf("ERROR: In xnl_integral. Memory allocation failed! Check variable weights.\n");
		exit(EXIT_FAILURE);
	}

    double total, left, right;

    /* Initialize roots and weights */

    create(roots, weights, points);

    /* Makes the subs x = tan(u) to avoid INFINITY intervals
       Remember: tan(+-inf) = +-pi/2
    */   
    if (src_isinf(a) || src_isinf(b))
    {
        if (src_isinf(a))
        {
            if (a > 0)
            {
                a = pi/2;
            } else 
            {
                a = -pi/2;
            }
        }
        else
        {
            a = atan(a); /* u = atan(x) */
        }

        if (src_isinf(b))
        {
            if (b > 0)
            {
                b = pi/2;
            } else 
            {
                b = -pi/2;
            }
        }
        else
        {
            b = atan(b); /* u = atan(x) */
        }

        func_pointer = f; /* Points to the passed function */

        if (!tolerance) 
        { 
       		/* nonadaptive approach */
       		total = GaussLegendre(g, a, b, points, roots, weights);
       		left = GaussLegendre(g, a, (a+b)/2, points, roots, weights);
       		right = GaussLegendre(g, (a+b)/2, b, points, roots, weights);

       		est_error = fabs(total - (left + right));

            retval = left + right;
        } 
        else /* adaptive approach */
        {
        	total = GaussLegendre(g, a, b, points, roots, weights);
            retval = Adaptive_Quad(g, a, b, points, total, tolerance, remaining_it, &est_error, roots, weights);
        }
    }
    else 
    {
        if (!tolerance) 
        { 
        	/* nonadaptive approach */
            total = GaussLegendre(f, a, b, points, roots, weights);
       		left = GaussLegendre(f, a, (a+b)/2, points, roots, weights);
       		right = GaussLegendre(f, (a+b)/2, b, points, roots, weights);

       		est_error = fabs(total - (left + right));

            retval = left + right;
        } 
        else /* adaptive approach */
        {
        	total = GaussLegendre(f, a, b, points, roots, weights);
            retval = Adaptive_Quad(f, a, b, points, total, tolerance, remaining_it, &est_error, roots, weights); 
        }
    }

    /* Free memory */

    free(roots);
    free(weights);

    /* Check if the quadrature worked */
    if (src_isnan(retval))
    {
        printf("Numerical integration failed. Possibly because the integral is divergent.\n");
    }

    if (error) *error = est_error;

    return retval;
}

/**********************************************
* Top level function
***********************************************/

static double integral(double (*f)(double), double a, double b, int points, double tolerance, int depth, double *error)
{
	/* Error handler */
	char error_msg[128] = "\0";
	int error_flag = 0;

	if (points <= 0)
	{
		strcat(error_msg,"ERROR: In integral. The number of points must be a natural number.\n");
		error_flag = 1;
	}

	if (tolerance < 0)
	{
		strcat(error_msg,"ERROR: In integral. The tolerance can't be a negative number.\n");
		error_flag = 1;
	}

	if (depth <= 0)
	{
		strcat(error_msg,"ERROR: In integral. The depth must be a positive integer.\n");
		error_flag = 1;
	}

	if (error_flag)
	{
		puts(error_msg);
		exit(EXIT_FAILURE);
	}

	return xnl_integral(f, a, b, points, tolerance, depth, error);
}

/**********************************************
* Python interface
***********************************************/

/*callback function*/

/*
double xnl_integral(double (*)(double),double,double,const unsigned,double,unsigned,double *);
*/

static PyObject *int_cb_callable;

/*Implementation of double (*func)(double)*/

static double int_callback(double callback)
{
	PyObject *retval;
	double result; /*returned value*/

	PyObject *value = Py_BuildValue("d", callback);
	PyObject *arglist = PyTuple_New(1);

	PyTuple_SetItem(arglist, 0, value); /*steals the reference of value*/

	/*call the python function/object saved below*/
	retval = PyObject_CallObject(int_cb_callable, arglist);

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
PyObject *py_integral(PyObject *self, PyObject *args, PyObject *kwargs)
{
	/*arguments*/
	double a;
	double b;
	unsigned points = 3;
	double tolerance = 1e-9;
	unsigned remaining_it = 10;
	/*values to be returned*/
	double error = 0;
	double result;
	PyObject *ret_tuple = PyTuple_New(2); 
	PyObject *ret_integral;
	PyObject *ret_error;

	static char *kwlist[] = {"f", "a", "b", "points", "tolerance", "depth", NULL};

	/* parse the input tuple with keywords */

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Odd|$IdI", kwlist, &int_cb_callable, &a, &b, &points, &tolerance, &remaining_it))
	{
		return NULL;
	}

	/* ensure the first parameter is a callable */
	if (!PyCallable_Check(int_cb_callable))
	{
		PyErr_SetString(PyExc_TypeError, "integral: a callable is required\n");
		return NULL;
	}

	/* Call the integral function */
	result = integral(int_callback,a,b,points,tolerance,remaining_it,&error);

	ret_integral = Py_BuildValue("d", result);
	ret_error = Py_BuildValue("d", error);

	PyTuple_SetItem(ret_tuple, 0, ret_integral); /*steals the reference of ret_integral*/
	PyTuple_SetItem(ret_tuple, 1, ret_error); /*steals the reference of ret_error*/

	return ret_tuple;
}
