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
	/* Francesco Tricomi formula */
	int k = (n / 2) + (n % 2);
	double dk = (double) k;
	double dn = (double) n;
	double theta = XNL_PI * (4 * dk - 1.0) / (4 * dn + 2.0);
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

	while (XNL_EPS < fabs(step) && l < 10)
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
		theta = XNL_PI * (4 * dk - 1.0) / (4 * dn + 2.0);
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

		while (XNL_EPS < fabs(step) && l < 10)
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

static double GaussLegendre(void (*f)(xparrayObject *, PyObject *, double *, int, double *, double, double), 
							xparrayObject *arg_array,
							PyObject *f_args, 
							double *vec, 
							double a, 
							double b, 
							int points, 
							double *roots, 
							double *weights)
{
	double c1 = (b - a)/2, c2 = (b + a)/2;
    double total = 0;
    int i;

    /*Vectorized function*/
    f(arg_array, f_args, vec, points, roots, c1, c2);

    for (i = 0; i < points; ++i) 
    {
        /*total += weights[i] * f(c1 * roots[i] + c2);*/
        total += weights[i] * vec[i];
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

static double Adaptive_Quad(void (*f)(xparrayObject *, PyObject *, double *, int, double *, double, double), /* function to integrate */
                   	 		xparrayObject *arg_array,
                   	 		PyObject *f_args, 
                   	 		double *vec, /*vector passed to the vectorized f function*/
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

	/* Checks if there is remaining iterations or if an error occurred*/
	if (!remaining_it || PyErr_Occurred())
	{
	    *error = XNL_NAN;

	    return XNL_NAN;
	}

	left = GaussLegendre(f, arg_array, f_args, vec, a, (a+b)/2, points, roots, weights); /* integrate over the left interval */
	right = GaussLegendre(f, arg_array, f_args, vec, (a+b)/2, b, points, roots, weights); /* integrate over the right interval */

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
	    retvalL = Adaptive_Quad(f, arg_array, f_args, vec, a, (a + b)/2, points, left, tolerance, remaining_it - 1, error, roots, weights);

	    if (src_isnan(retvalL))
	    {
	        /* Checks if the left node isn't ok */
	        *error = XNL_NAN;
	    
	    	return retvalL;
	    } 
	    else 
	    {
	    	/* The left node is ok. Then the right node can be computed */
	        retvalR = Adaptive_Quad(f, arg_array, f_args, vec, (a + b)/2, b, points, right, tolerance, remaining_it - 1, error, roots, weights);

	    	return retvalL + retvalR; /* returns the left + right nodes */
	    }
	}
}

/**********************************************
* Adaptive gaussian quadrature
***********************************************/

static void (*func_pointer)(xparrayObject *, PyObject *, double *, int, double *, double, double);

static void g(xparrayObject *arg_array, PyObject *f_args, double *vec, int points, double *roots, double c1, double c2)
{
	double *new_x = malloc(points * sizeof(double));

	if (new_x == NULL)
	{
		printf("ERROR: In g. Memory allocation failed. Please check variable new_x.\n");
		exit(EXIT_FAILURE);
	}

	int i;
	for (i = 0; i < points; i++)
	{
		new_x[i] = tan(c1*roots[i] + c2);
	}

	func_pointer(arg_array, f_args, vec, points, new_x, 1, 0); 

	double cos_new_x;
	for (i = 0; i < points; i++)
	{
		cos_new_x = cos(c1 * roots[i] + c2);
		vec[i] = vec[i] / (cos_new_x * cos_new_x);
	}

	free(new_x);
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

static double xnl_integral(void (*f)(xparrayObject *, PyObject *, double *, int, double *, double, double), /* function to integrate */
               	 		   PyObject *f_args, 
               	 		   double a, /* left interval boundary */
               			   double b, /* right interval boundary */
						   int points, /* number of points for interpolation */
               			   double tolerance, /* error tolerance */
               			   int remaining_it, /* maximum number of iterations per division */
               			   double* error) /* stores the error value */
{
	xparrayObject *arg_array = NULL;

	arg_array = PyXParray_New(1, points);

	/* Set the first element of the tuple f_args to be the array argument (x) */
	/* Insert a referente to arg_array at the position 0 */
	if (f_args)
		PyTuple_SetItem(f_args, 0, (PyObject *) arg_array); /* steals the reference of arg_array */

    double retval; /* Control variable */

    double est_error = 0; /* Estimated error */

    /* Legendre Polynomial roots and weights */

    double *roots = malloc(points * sizeof(double));
    double *weights = malloc(points * sizeof(double));
    double *vec = malloc(points * sizeof(double)); /*to vectorize f function*/

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

	if (vec == NULL)
	{
		printf("ERROR: In xnl_integral. Memory allocation failed! Check variable vec.\n");
		exit(EXIT_FAILURE);
	}	

    double total, left, right;

    void (*func)(xparrayObject *, PyObject *, double *, int, double *, double, double);

    func = f;

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
                a = XNL_PI/2;
            } else 
            {
                a = -XNL_PI/2;
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
                b = XNL_PI/2;
            } else 
            {
                b = -XNL_PI/2;
            }
        }
        else
        {
            b = atan(b); /* u = atan(x) */
        }

        func_pointer = f; /* Points to the passed function */
        func = g;
    }
    
    if (!tolerance) 
    { 
    	/* nonadaptive approach */
        total = GaussLegendre(func, arg_array, f_args, vec, a, b, points, roots, weights);
   		left = GaussLegendre(func, arg_array, f_args, vec, a, (a+b)/2, points, roots, weights);
   		right = GaussLegendre(func, arg_array, f_args, vec, (a+b)/2, b, points, roots, weights);

   		est_error = fabs(total - (left + right));

        retval = left + right;
    } 
    else /* adaptive approach */
    {
    	total = GaussLegendre(func, arg_array, f_args, vec, a, b, points, roots, weights);
        retval = Adaptive_Quad(func, arg_array, f_args, vec, a, b, points, total, tolerance, remaining_it, &est_error, roots, weights); 
    }

    /* Free memory */

    free(roots);
    free(weights);
    free(vec);

    /* Check if the quadrature worked */
    if (src_isnan(retval) && !PyErr_Occurred())
    {
        printf("WARNING: In xnl_integral. Numerical integration failed. Possibly because the integral is divergent.\n");
    }

    if (error) *error = est_error;

	if (!f_args)
		Py_XDECREF(arg_array);

    return retval;
}

/**********************************************
* Top level function
***********************************************/

static double integral(void (*f)(xparrayObject *, PyObject *, double *, int, double *, double, double), 
					   PyObject *f_args, 
					   double a, 
					   double b, 
					   int points, 
					   double tolerance, 
					   int depth, 
					   double *error)
{
	/* Error handler */
	char error_msg[128] = "\0";
	int error_flag = 0;

	/* Handles error of args parameter */	
	if (PyErr_Occurred())
	{
		return 0; /* just returning anything */
	}

	if (points <= 0)
	{
		strcat(error_msg,"The number of points must be a natural number. ");
		error_flag = 1;
	}

	if (tolerance < 0)
	{
		strcat(error_msg,"The tolerance can't be a negative number. ");
		error_flag = 1;
	}

	if (depth <= 0)
	{
		strcat(error_msg,"The depth must be a positive integer. ");
		error_flag = 1;
	}

	if (error_flag)
	{
		PyErr_SetString(PyExc_TypeError, error_msg); /* This will set error flag */
		return 0; /* just returning anything */
	}

	return xnl_integral(f, f_args, a, b, points, tolerance, depth, error);
}

/**********************************************
* Python interface
***********************************************/

/*callback function*/

/*
double integral(void (*)(xparrayObject *, PyObject *, double *, int, double *, double, double),
				PyObject *, 
				double,
				double,
				int, double, int,
				double *);
*/

static PyObject *int_cb_callable = NULL;

/*Implementation of void (*func)(xparrayObject *, PyObject *, double *, int, double *, double, double)*/

static void int_callback(xparrayObject *arg_array,
						 PyObject *f_args, 
						 double *vec, 
						 int points, 
						 double *roots, 
						 double c1, 
						 double c2)
{
	unsigned i;
    for (i = 0; i < points; i++)
	{
		arg_array->data[0][i] = c1 * roots[i] + c2;
	}

    PyObject *ret_val;
	/*Returns a reference to ret_val*/
	if (!f_args)
	{
		ret_val = PyObject_CallFunctionObjArgs(int_cb_callable, (PyObject *) arg_array, NULL);
	}
	else 
	{
		/* user_data available */
		ret_val = PyObject_CallObject(int_cb_callable, f_args);
	}

	if (ret_val && !PyXParray_Check(ret_val))
	{
	   	xparrayObject *ret_array = (xparrayObject *) ret_val;
		/*vec = ret_array->data[0];*/ /*vec now points to ret_array->data[0]*/
		for (i = 0; i < points; i++)
		{
			vec[i] = ret_array->data[0][i];
		}
	}
	else 
	{
		/* this will set the error flag */
		PyErr_SetString(PyExc_RuntimeError, "Callable object int: bad behavior");
	}

	Py_XDECREF(ret_val);
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
	PyObject *fc = PyObject_GetAttrString(int_cb_callable, "__code__");
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
						"integral: args not available for this integrand. "
						"It must have more than one argument.");
		return NULL;
	}

	if (count > 1 && !user_data)
	{
		PyErr_SetString(PyExc_TypeError, 
						"integral: invalid number of arguments. The integrand "
						"appears to have more than one argument, so you must "
						"use \"args\" to pass the other arguments to it.");
		return NULL;
	}

	if (!user_data)
		return NULL;

	/* ensure the parameter args is a tuple */
	if (!PyTuple_Check(user_data))
	{
		PyErr_SetString(PyExc_TypeError, "integral: args must be a tuple");
		return NULL;
	}

	/* ensure args is a tuple of numbers */
	Py_ssize_t tuplen = PyTuple_Size(user_data);

	if (tuplen != (count - 1))
	{
		PyErr_SetString(PyExc_RuntimeError, "integral: args' size must be equal to "
											"the number of extra arguments of the "
											"integrand.");
		return NULL;
	}

	Py_ssize_t i;

	for (i = 0; i < tuplen; i++)
	{
		if (!PyNumber_Check(PyTuple_GetItem(user_data, i)))
		{
			PyErr_SetString(PyExc_TypeError, "integral: args must be a tuple of numbers");
			return NULL;
		}
	}

	/* user_data is ok. Now creates a tuple with the callable args */
	/**
	 * The first argument of f is a xparray. So it's size is 1 + tuplen.
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

/*****************************************************************************
 * 
 * 							Top level function
 * 
 *****************************************************************************
 */
PyObject *PyXP_Integral(PyObject *self, PyObject *args, PyObject *kwargs)
{
	/*integral arguments*/
	double a;
	double b;
	PyObject *user_data = NULL;
	int points = 3;
	double tolerance = 1e-9;
	int remaining_it = 10;
	/*values to be returned*/
	double error = 0;
	double result;

	static char *kwlist[] = {"f", "a", "b", "args", "points", "tolerance", "depth", NULL};

	/* parse the input tuple with keywords */

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Odd|$OIdI", kwlist, &int_cb_callable, &a, &b, &user_data, &points, &tolerance, &remaining_it))
	{
		PyErr_SetString(PyExc_TypeError, "integral: failed to parse arguments");
		return NULL;
	}

	/* ensure the first parameter is a callable */
	if (!PyCallable_Check(int_cb_callable))
	{
		PyErr_SetString(PyExc_TypeError, "integral: a callable is required");
		return NULL;
	}

	PyObject *f_args = check_args(user_data);

	/* Calls the integral function */
	result = integral(int_callback,f_args,a,b,points,tolerance,remaining_it,&error);

	Py_XDECREF(f_args);

	/* Error occurred when working with integral.
	 * The only possibility is error due to a bad 
	 * behavior of the callable object. Just returns 
	 * NULL to indicate to python that an error occurred 
	 * without closing the interpreter.
	 */ 
	if (PyErr_Occurred())
	{
		return NULL;
	}

	/* Everything ok, then returns the tuple (result, error) */

	PyObject *ret_tuple = PyTuple_New(2); 
	PyObject *ret_integral;
	PyObject *ret_error;

	ret_integral = Py_BuildValue("d", result);
	ret_error = Py_BuildValue("d", error);

	PyTuple_SetItem(ret_tuple, 0, ret_integral); /*steals the reference of ret_integral*/
	PyTuple_SetItem(ret_tuple, 1, ret_error); /*steals the reference of ret_error*/

	return ret_tuple;
}
