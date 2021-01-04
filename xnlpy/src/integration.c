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

static double GaussLegendre(void (*f)(xparrayObject *, double *, int, double *, double, double), 
							xparrayObject *arg_array,
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
    f(arg_array, vec, points, roots, c1, c2);

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

static double Adaptive_Quad(void (*f)(xparrayObject *, double *, int, double *, double, double), /* function to integrate */
                   	 		xparrayObject *arg_array,
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

	/* Checks if there is remaining iterations */
	if (!remaining_it)
	{
	    *error = XNL_NAN;

	    return XNL_NAN;
	}

	left = GaussLegendre(f, arg_array, vec, a, (a+b)/2, points, roots, weights); /* integrate over the left interval */
	right = GaussLegendre(f, arg_array, vec, (a+b)/2, b, points, roots, weights); /* integrate over the right interval */

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
	    retvalL = Adaptive_Quad(f, arg_array, vec, a, (a + b)/2, points, left, tolerance, remaining_it - 1, error, roots, weights);

	    if (src_isnan(retvalL))
	    {
	        /* Checks if the left node isn't ok */
	        *error = XNL_NAN;
	    
	    	return retvalL;
	    } 
	    else 
	    {
	    	/* The left node is ok. Then the right node can be computed */
	        retvalR = Adaptive_Quad(f, arg_array, vec, (a + b)/2, b, points, right, tolerance, remaining_it - 1, error, roots, weights);

	    	return retvalL + retvalR; /* returns the left + right nodes */
	    }
	}
}

/**********************************************
* Adaptive gaussian quadrature
***********************************************/

static void (*func_pointer)(xparrayObject *, double *, int, double *, double, double);

static void g(xparrayObject *arg_array, double *vec, int points, double *roots, double c1, double c2)
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

	func_pointer(arg_array, vec, points, new_x, 1, 0); 

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

static double xnl_integral(void (*f)(xparrayObject *, double *, int, double *, double, double), /* function to integrate */
               	 		   double a, /* left interval boundary */
               			   double b, /* right interval boundary */
						   int points, /* number of points for interpolation */
               			   double tolerance, /* error tolerance */
               			   int remaining_it, /* maximum number of iterations per division */
               			   double* error) /* stores the error value */
{
	xparrayObject *arg_array = NULL;

	arg_array = PyXParray_New(1, points);

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

    void (*func)(xparrayObject *, double *, int, double *, double, double);

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
        total = GaussLegendre(func, arg_array, vec, a, b, points, roots, weights);
   		left = GaussLegendre(func, arg_array, vec, a, (a+b)/2, points, roots, weights);
   		right = GaussLegendre(func, arg_array, vec, (a+b)/2, b, points, roots, weights);

   		est_error = fabs(total - (left + right));

        retval = left + right;
    } 
    else /* adaptive approach */
    {
    	total = GaussLegendre(func, arg_array, vec, a, b, points, roots, weights);
        retval = Adaptive_Quad(func, arg_array, vec, a, b, points, total, tolerance, remaining_it, &est_error, roots, weights); 
    }
    
    /* Free memory */

    free(roots);
    free(weights);
    free(vec);

    /* Check if the quadrature worked */
    if (src_isnan(retval))
    {
        printf("ERROR: In xnl_integral. Numerical integration failed. Possibly because the integral is divergent.\n");
    }

    if (error) *error = est_error;

	Py_XDECREF(arg_array);

    return retval;
}

/**********************************************
* Top level function
***********************************************/

static double integral(void (*f)(xparrayObject *, double *, int, double *, double, double), 
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

/**
 * Double integral numerical approximation
 *
 * TwoD algorithm by L. F. Shampine
 */

double c_scalar = 0;
double d_scalar = 0;
int c_callable = 0;
int d_callable = 0;

static void phivar(PyObject *f, xparrayObject *x, double *dst, int flag, double scalar)
{
	unsigned i;

	if (!flag)
	{
		for (i = 0; i < 14; i++)
		{
			dst[i] = scalar;
		}

		return;
	}

	PyObject *ret_val;
	/*Returns a reference to ret_val*/
	ret_val = PyObject_CallFunctionObjArgs(f, (PyObject *) x, NULL);

	if (ret_val && !PyXParray_Check(ret_val))
	{
	   	xparrayObject *ret_array = (xparrayObject *) ret_val;
		/*vec = ret_array->data[0];*/ /*vec now points to ret_array->data[0]*/
		for (i = 0; i < 14; i++)
		{
			dst[i] = ret_array->data[0][i];
		}
	}
	else 
	{
		PyErr_SetString(PyExc_TypeError, "Callable object phivar: bad behavior");
	}

	Py_XDECREF(ret_val);
}

static void tensor( double *Qsub, double *esub, 
					double thetaL, double thetaR, double phiB, double phiT,
					double (*NArray)[2],
					double *wt3, double *wt7,
					void (*func)(xparrayObject *, xparrayObject *, xparrayObject *),
					xparrayObject *z_array, 
					xparrayObject *x_array, 
					xparrayObject *y_array,
					double a, double b,
					PyObject *c, PyObject *d,
					xparrayObject *x,
					int Singular)
{
	double dtheta = thetaR - thetaL;

	/**
	 * theta = bsxfun(@plus,dtheta*NARRAY,thetaL+dtheta*[0.25,0.75]);
     * theta = reshape(theta,1,[]); 
	 */
	double theta[14];

/*
	int i;

	for (i = 0; i < 14; i++)
	{
		theta[i] = (i < 7) ? dtheta * (NArray[i][0] + 0.25) + thetaL : dtheta * (NArray[i-7][1] + 0.75) + thetaL;
	}
*/

/*
	double x[14];

	for (i = 0; i < 14; i++)
	{
		x[i] = (Singular) ? 0.5*(b + a) + 0.5*(b - a)*cos(theta[i]) : theta[i]; 
	}
*/

/*
	double xmat[14][14];

	for (i = 0; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
		{
			xmat[i][j] = x[j];
		}
	}
*/
	double dphi = phiT - phiB;

	/**
	 * phi = bsxfun(@plus,dphi*NARRAY,phiB+dphi*[0.25,0.75]);
     * phi = reshape(phi,1,[])';  
	 */
	
	double phi[14];
/*
	for (i = 0; i < 14; i++)
	{
		phi[i] = (i < 7) ? dphi * (NArray[i][0] + 0.25) + phiB : dphi * (NArray[i-7][1] + 0.75) + phiB;
	}
*/
	double top[14], bottom[14], dydt[14];
/*
	for (i = 0; i < 14; i++)
	{
		top[i] = phivar(d, x[i]);
		bottom[i] = phivar(c, x[i]);
		dydt[i] = top[i] - bottom[i];
	}
*/
	double t[14];
/*
	for (i = 0; i < 14; i++)
	{
		t[i] = (Singular) ? 0.5 * (1 + cos(phi[i])) : phi[i];
	}
*/

/*
	double ymat[14][14];

	for (i = 0; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
		{
			ymat[i][j] = bottom[j] + t[i] * dydt[j];
		}
	}
*/

	double temp[14][14];
/*
	for (i = 0; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
		{
			temp[i][j] = (Singular) ? 0.25 * (b - a) * sin(phi[i]) * (dydt[j] * sin(theta[j])) : dydt[j];
		}
	}
*/

/*
	double zmat[14][14];

	for (i = 0; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
		{
			zmat[i][j] = func(xmat[i][j], ymat[i][j]) * temp[i][j];
		}
	}
*/

	unsigned i, j;

	for (i = 0; i < 14; i++)
	{
		theta[i] = (i < 7) ? dtheta * (NArray[i][0] + 0.25) + thetaL : dtheta * (NArray[i-7][1] + 0.75) + thetaL;
		x->data[0][i] = (Singular) ? 0.5*(b + a) + 0.5*(b - a)*cos(theta[i]) : theta[i]; 
		x_array->data[0][i] = x->data[0][i];
		phi[i] = (i < 7) ? dphi * (NArray[i][0] + 0.25) + phiB : dphi * (NArray[i-7][1] + 0.75) + phiB;
		t[i] = (Singular) ? 0.5 * (1 + cos(phi[i])) : phi[i];
	}

	phivar(c, x, bottom, c_callable, c_scalar);
	phivar(d, x, top, d_callable, d_scalar);

	for (i = 0; i < 14; i++)
	{
		/*
		theta[i] = (i < 7) ? dtheta * (NArray[i][0] + 0.25) + thetaL : dtheta * (NArray[i-7][1] + 0.75) + thetaL;
		x[i] = (Singular) ? 0.5*(b + a) + 0.5*(b - a)*cos(theta[i]) : theta[i]; 
		x_array->data[0][i] = x[i];
		phi[i] = (i < 7) ? dphi * (NArray[i][0] + 0.25) + phiB : dphi * (NArray[i-7][1] + 0.75) + phiB;
		top[i] = phivar(d, x[i]);
		bottom[i] = phivar(c, x[i]);
		*/
		dydt[i] = top[i] - bottom[i];
		/*t[i] = (Singular) ? 0.5 * (1 + cos(phi[i])) : phi[i];*/
		y_array->data[0][i] = bottom[i] + t[0] * dydt[i];
		temp[0][i] = (Singular) ? 0.25 * (b - a) * sin(phi[0]) * (dydt[i] * sin(theta[i])) : dydt[i];
		/*zmat[0][i] = func(xmat[0][i], ymat[0][i]) * temp[0][i];*/
	}

	for (i = 1; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
		{
			x_array->data[i][j] = x->data[0][j];
			y_array->data[i][j] = bottom[j] + t[i] * dydt[j];
			temp[i][j] = (Singular) ? 0.25 * (b - a) * sin(phi[i]) * (dydt[j] * sin(theta[j])) : dydt[j];
			/*zmat[i][j] = func(xmat[i][j], ymat[i][j]) * temp[i][j];*/
		}
	}

	func(z_array, x_array, y_array);

	for (i = 0; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
			z_array->data[i][j] *= temp[i][j];
	}

	/* reset Qsub and esub */
	for (i = 0; i < 4; i++)
	{
		Qsub[i] = esub[i] = 0;
	}

	double Jdthetadphi = (dtheta * dphi) / 16;

	/* Tensor product Gauss 3 point formula */
	for (i = 0; i < 7; i++)
	{
		double soma0 = 0,
			   soma1 = 0,
			   soma2 = 0,
			   soma3 = 0;
		for (j = 0; j < 7; j++)
		{
			soma0 += wt3[j] * z_array->data[j][i];
			soma1 += wt3[j] * z_array->data[j][i+7];
			soma2 += wt3[j] * z_array->data[j+7][i];
			soma3 += wt3[j] * z_array->data[j+7][i+7];
		}
		esub[0] += wt3[i] * soma0 * Jdthetadphi;
		esub[1] += wt3[i] * soma1 * Jdthetadphi;
		esub[2] += wt3[i] * soma2 * Jdthetadphi;
		esub[3] += wt3[i] * soma3 * Jdthetadphi;
	}

	/* Tensor product Kronrod 7 point formula */
	for (i = 0; i < 7; i++)
	{
		double soma0 = 0,
			   soma1 = 0,
			   soma2 = 0,
			   soma3 = 0;
		for (j = 0; j < 7; j++)
		{
			soma0 += wt7[j] * z_array->data[j][i];
			soma1 += wt7[j] * z_array->data[j][i+7];
			soma2 += wt7[j] * z_array->data[j+7][i];
			soma3 += wt7[j] * z_array->data[j+7][i+7];
		}
		Qsub[0] += wt7[i] * soma0 * Jdthetadphi;
		Qsub[1] += wt7[i] * soma1 * Jdthetadphi;
		Qsub[2] += wt7[i] * soma2 * Jdthetadphi;
		Qsub[3] += wt7[i] * soma3 * Jdthetadphi;
	}

	for (i = 0; i < 4; i++)
	{
		esub[i] = fabs(esub[i] - Qsub[i]);
	}
}

static double Save2List(double *Qsub, double *esub,
				 		double thetaL, double thetaR, double phiB, double phiT,
			  	 		double ***LList, int *NList,
			  	 		double tol, double *err_ok, double adjust, double Area,
			  	 		unsigned *sizelist)
{
	double dtheta = thetaR - thetaL;
	double dphi = phiT - phiB;
	double localtol = tol * (dtheta/2) * (dphi/2) / Area;
	double Q = 0;
	double **List = *LList;

	int i;
	for (i = 0; i < 4; i++)
		Q += Qsub[i];

	double tolval = 100*XNL_EPS*fabs(Q);

	/*Não poderia usar tol?*/
	/*Provalmente Qsub e esub são constantemente atualizados*/

	localtol = (localtol > tolval) ? localtol : tolval;

	double adjerr[4];

	for (i = 0; i < 4; i++)
	{
		adjerr[i] = adjust * esub[i];
	} 

	int j;
	if (*NList + 4 > *sizelist)
	{
		unsigned old_sizelist = *sizelist;
		/*Allocate more memory to List*/
		*sizelist += 100;
		double **temp = realloc(List, *sizelist * sizeof(double *));
		if (temp == NULL)
		{
			printf("ERROR: In TwoD. Memory allocation failed after realloc. Please check the variable temp.\n");
			exit(EXIT_FAILURE);
		}

		List = temp;
		for (i = old_sizelist; i < *sizelist; i++)
		{
			List[i] = malloc(7 * sizeof(double));
			if (List[i] == NULL)
			{
				printf("ERROR: In TwoD. Memory allocation failed after realloc. Please check the variable List[i].\n");
				exit(EXIT_FAILURE);
			}
			for (j = 0; j < 7; j++)
			{
				List[i][j] = 0;
			}
		}
	}

	if (adjerr[0] > localtol)
	{
		*NList += 1;
		
		List[*NList-1][0] = Qsub[0];
		List[*NList-1][1] = esub[0];
		List[*NList-1][2] = thetaL;
		List[*NList-1][3] = thetaL + dtheta/2;
		List[*NList-1][4] = phiB;
		List[*NList-1][5] = phiB + dphi/2;
		List[*NList-1][6] = adjerr[0];
	}
	else
	{
		*err_ok += adjerr[0];
	}

	if (adjerr[1] > localtol)
	{
		*NList += 1;
		
		List[*NList-1][0] = Qsub[1];
		List[*NList-1][1] = esub[1];
		List[*NList-1][2] = thetaL + dtheta/2;
		List[*NList-1][3] = thetaR;
		List[*NList-1][4] = phiB;
		List[*NList-1][5] = phiB + dphi/2;
		List[*NList-1][6] = adjerr[1];
	}
	else
	{
		*err_ok += adjerr[1];
	}

	if (adjerr[2] > localtol)
	{
		*NList += 1;
		
		List[*NList-1][0] = Qsub[2];
		List[*NList-1][1] = esub[2];
		List[*NList-1][2] = thetaL;
		List[*NList-1][3] = thetaL + dtheta/2;
		List[*NList-1][4] = phiB + dphi/2;
		List[*NList-1][5] = phiT;
		List[*NList-1][6] = adjerr[2];
	}
	else
	{
		*err_ok += adjerr[2];
	}

	if (adjerr[3] > localtol)
	{
		*NList += 1;

		List[*NList-1][0] = Qsub[3];
		List[*NList-1][1] = esub[3];
		List[*NList-1][2] = thetaL + dtheta/2;
		List[*NList-1][3] = thetaR;
		List[*NList-1][4] = phiB + dphi/2;
		List[*NList-1][5] = phiT;
		List[*NList-1][6] = adjerr[3];
	}
	else
	{
		*err_ok += adjerr[3];
	}
	
	double errbnd = 0;

	for (i = 0; i < *sizelist; i++)
	{
		errbnd += List[i][6];
	}

	*LList = List;

	return errbnd + *err_ok;
}

static void NextEntry(double *q, double *e, 
			   		  double *thetaL, double *thetaR,
			   		  double *phiB, double *phiT,
			   		  double **List, int *NList,
			   		  unsigned *sizelist)
{
	double max = fabs(List[0][6]);
	unsigned max_index = 0;

	unsigned i;
	for (i = 1; i < *NList; i++)
	{
		double absList = fabs(List[i][6]);
		if (absList > max)
		{
			max = absList;
			max_index = i;
		}
	}

	*q = List[max_index][0];
	*e = List[max_index][1];
	*thetaL = List[max_index][2];
	*thetaR = List[max_index][3];
	*phiB = List[max_index][4];
	*phiT = List[max_index][5];

	*NList = *NList - 1;

	unsigned j;
	for (i = max_index; i < *sizelist-1; i++)
	{
		for (j = 0; j < 7; j++)
		{
			List[i][j] = List[i+1][j];
		}
	}

	free(List[*sizelist-1]);

	*sizelist = *sizelist - 1;
}

static double max(double x, double y)
{
	return (x > y) ? x : y;
}

static void (*f_ptr)(xparrayObject *, xparrayObject *, xparrayObject *);
static void (*f2_ptr)(xparrayObject *, xparrayObject *, xparrayObject *);
static xparrayObject *x_aux = NULL;
static xparrayObject *y_aux = NULL;

static void polar(xparrayObject *ret_array, xparrayObject *theta, xparrayObject *r)
{
	unsigned i, j;

	for (i = 0; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
		{
			x_aux->data[i][j] = r->data[i][j] * cos(theta->data[i][j]);
			y_aux->data[i][j] = r->data[i][j] * sin(theta->data[i][j]);
		}
	}

	/* f_ptr(r * cos(theta), r * sin(theta)) * r */
	f_ptr(ret_array, x_aux, y_aux);

	for (i = 0; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
		{
			ret_array->data[i][j] *= r->data[i][j];
		}
	}
}

static void func_inf(xparrayObject *ret_array, xparrayObject *x, xparrayObject *y)
{
	unsigned i, j;

	for (i = 0; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
		{
			x_aux->data[i][j] = tan(x->data[i][j]);
			y_aux->data[i][j] = tan(y->data[i][j]);
		}
	}

	f2_ptr(ret_array, x_aux, y_aux);

	for (i = 0; i < 14; i++)
	{
		for (j = 0; j < 14; j++)
		{
			ret_array->data[i][j] /= (cos(x->data[i][j]) * cos(x->data[i][j]) * cos(y->data[i][j]) * cos(y->data[i][j]));
		}
	}
}

static double TwoD(void (*f)(xparrayObject *, xparrayObject *, xparrayObject *), 
				   double a, double b,
				   PyObject *c, PyObject *d,
				   double AbsTol, /*1e-5*/
				   double RelTol, /*0*/
				   int Sector, /*false*/
				   int Singular, /*false*/
				   double *errbnd)
{
	xparrayObject *x = NULL;
	xparrayObject *x_array = NULL;
	xparrayObject *y_array = NULL;
	xparrayObject *z_array = NULL;

	x = PyXParray_New(1,14);
	x_array = PyXParray_New(14, 14);
	y_array = PyXParray_New(14, 14);
	z_array = PyXParray_New(14, 14);
	x_aux = PyXParray_New(14, 14);
	y_aux = PyXParray_New(14, 14);

	unsigned sizelist = 200;

	void (*func)(xparrayObject *, xparrayObject *, xparrayObject *);
	f_ptr = f; /* Points to the passed function */

	if (Sector)
	{
		func = polar;	
	}
	else
	{
		func = f;
	}

	/* Makes the subs x = tan(u) to avoid INFINITY intervals
       Remember: tan(+-inf) = +-pi/2
    */   
    if (src_isinf(a) || src_isinf(b) || src_isinf(c_scalar) || src_isinf(d_scalar))
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

        if (src_isinf(c_scalar))
        {
            if (c_scalar > 0)
            {
                c_scalar = XNL_PI/2;
            } else 
            {
                c_scalar = -XNL_PI/2;
            }
        }
        else
        {
            c_scalar = atan(c_scalar); /* u = atan(x) */
        }

        if (src_isinf(d_scalar))
        {
            if (d_scalar > 0)
            {
                d_scalar = XNL_PI/2;
            } else 
            {
                d_scalar = -XNL_PI/2;
            }
        }
        else
        {
            d_scalar = atan(d_scalar); /* u = atan(x) */
        }

        f2_ptr = (Sector) ? polar : f;;
        func = func_inf;
    }

	double thetaL = 0,
		   thetaR = XNL_PI,
		   phiB = 0,
		   phiT = XNL_PI;
	if (!Singular)
	{
		thetaL = a;
		thetaR = b;
		phiB = 0;
		phiT = 1;
	}

	double Area = (thetaR - thetaL) * (phiT - phiB);

	/*
	double Nodes[7] = {-0.9604912687080202, -0.7745966692414834, -0.4342437493468026,
             			0, 0.4342437493468026,  0.7745966692414834,  0.9604912687080202};
    */
    /*int NNodes = 7;*/
    /*int OneVec[14] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};*/
    /*NArray = 0.25 * [Nodes(:), Nodes(:)]*/
    double NArray[7][2] = {{-0.240122817177005,  -0.240122817177005},
						   {-0.193649167310371,  -0.193649167310371},
						   {-0.108560937336701,  -0.108560937336701},
						   {                 0,                   0},
						   { 0.108560937336701,   0.108560937336701},
						   { 0.193649167310371,   0.193649167310371},
						   { 0.240122817177005,   0.240122817177005}};

	double wt3[7] = {0, 5.0/9, 0, 8.0/9, 0, 5.0/9, 0};
	double wt7[7] = {0.1046562260264672, 0.2684880898683334, 0.4013974147759622,
		             0.4509165386584744, 0.4013974147759622, 0.2684880898683334,
		             0.1046562260264672};

	double Qsub[4] = {0},
		   esub[4] = {0};

	tensor(Qsub, esub, 
		   thetaL, thetaR, phiB, phiT, 
		   NArray, wt3, wt7, 
		   func, z_array, x_array, y_array, 
		   a, b, c, d, x, Singular);

	double Q = 0;

	int i;
	for (i = 0; i < 4; i++)
		Q += Qsub[i];

	double tol = 100 * XNL_EPS * fabs(Q);
	double err_ok = 0;
	int NList = 0;
	double adjust = 1;
	double **List;

	List = malloc(200 * sizeof(double *));
	if (List == NULL)
	{
		printf("ERROR: In TwoD. Memory allocation failed. Please check the variable List.\n");
		exit(EXIT_FAILURE);
	}
	int j;
	for (i = 0; i < 200; i++)
	{
		List[i] = malloc(7 * sizeof(double));
		if (List[i] == NULL)
		{
			printf("ERROR: In TwoD. Memory allocation failed. Please check the variable List[i].\n");
			exit(EXIT_FAILURE);
		}
		for (j = 0; j < 7; j++)
		{
			List[i][j] = 0;
		}
	}

	*errbnd = Save2List(Qsub, esub, 
						thetaL, thetaR, phiB, phiT,
						&List, &NList,
						tol, &err_ok, adjust, Area, 
						&sizelist);

	if (NList == 0 || *errbnd <= tol)
	{
		return Q;
	}

	while (1)
	{
		double q, e;

		NextEntry(&q, &e, &thetaL, &thetaR, &phiB, &phiT, List, &NList, &sizelist);

		/*printf("%.15f\n%.15f\n%.15f\n%.15f\n%.15f\n%.15f\n", q, e, thetaL, thetaR, phiB, phiT);*/

		tensor(Qsub, esub, 
		   thetaL, thetaR, phiB, phiT, 
		   NArray, wt3, wt7, 
		   func, z_array, x_array, y_array,  
		   a, b, c, d, x, Singular);

		double Newq = 0;

		for (i = 0; i < 4; i++)
			Newq += Qsub[i];

		adjust = fabs((q - Newq)/e);
		adjust = (adjust < 1) ? adjust : 1;

		Q = Q + (Newq - q);

		tol = RelTol * fabs(Q);
		tol = ( (AbsTol > tol) ? AbsTol : tol ) / 8;
		double temp = 100 * XNL_EPS * fabs(Q);
		tol = (tol > temp) ? tol : temp;

		*errbnd = Save2List(Qsub, esub, 
						    thetaL, thetaR, phiB, phiT,
						    &List, &NList,
						    tol, &err_ok, adjust, Area, 
						    &sizelist);

		/* Test for convergence and failure */
		if (NList == 0 || *errbnd <= tol)
			break;

		if (NList >= 2000)
		{
			if (*errbnd > max(AbsTol, max(100 * XNL_EPS, RelTol) * fabs(Q)))
			{
				printf("ERROR: In TwoD. Maximum number of subintervals: Q does NOT pass error test.\n");
				exit(EXIT_FAILURE);
			}
			else
			{
				printf("*******************************************************************************\n"
					   "WARNING: In TwoD. Maximum number of subintervals: Q appears to pass error test.\n"
					   "*******************************************************************************\n");
			}
			break;
		}
	}

	for (i = 0; i < sizelist; i++)
	{
		free(List[i]);
	}

	free(List);

	Py_XDECREF(x);
	Py_XDECREF(x_array);
	Py_XDECREF(y_array);
	Py_XDECREF(z_array);

	Py_XDECREF(x_aux);
	Py_XDECREF(y_aux);	

	return Q;
}

/**********************************************
* Top level function
***********************************************/

static double integral2(void (*f)(xparrayObject *, xparrayObject *, xparrayObject *), 
					    double a, double b,
					    PyObject *c, PyObject *d,
					    double AbsTol, /*1e-5*/
					    double RelTol, /*0*/
					    int Sector, /*false*/
					    int Singular, /*false*/
					    double *errbnd)
{
	/* Error handler */
	char error_msg[128] = "\0";
	int error_flag = 0;

	if (AbsTol < 0)
	{
		strcat(error_msg,"ERROR: In integral2. The parameter AbsTol cannot be negative.\n");
		error_flag = 1;
	}

	if (RelTol < 0)
	{
		strcat(error_msg,"ERROR: In integral2. The parameter RelTol cannot be negative.\n");
		error_flag = 1;
	}

	if (error_flag)
	{
		puts(error_msg);
		exit(EXIT_FAILURE);
	}

	return TwoD(f, a, b, c, d, AbsTol, RelTol, Sector, Singular, errbnd);
}

/**********************************************
* Python interface
***********************************************/

/*callback function*/

/*
double integral(void (*)(xparrayObject *, double *, int, double *, double, double),
				double,
				double,
				int, double, int,
				double *);
*/

static PyObject *int_cb_callable = NULL;

/*Implementation of void (*func)(xparrayObject *, double *, int, double *, double, double)*/

static void int_callback(xparrayObject *arg_array,
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
	ret_val = PyObject_CallFunctionObjArgs(int_cb_callable, (PyObject *) arg_array, NULL);

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
		PyErr_SetString(PyExc_TypeError, "Callable object int: bad behavior");
	}

	Py_XDECREF(ret_val);
}

/* Python interface of integral */
static PyObject *PYintegral1(PyObject *self, PyObject *args, PyObject *kwargs)
{
	/*integral arguments*/
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

	/* Calls the integral function */
	result = integral(int_callback,a,b,points,tolerance,remaining_it,&error);

	ret_integral = Py_BuildValue("d", result);
	ret_error = Py_BuildValue("d", error);

	PyTuple_SetItem(ret_tuple, 0, ret_integral); /*steals the reference of ret_integral*/
	PyTuple_SetItem(ret_tuple, 1, ret_error); /*steals the reference of ret_error*/

	return ret_tuple;
}

/*
double integral2(void (*)(xparrayObject *, double *, int, double *, double, double),
				 double,
				 double,
				 void (*)(xparrayObject *),
				 void (*)(xparrayObject *),
				 double, double, int, int, 
				 double *);
*/

/* Implementation of void (*func)(xparrayObject *, xparrayObject *, xparrayObject *) */

static void int2_callback(xparrayObject *z_array, xparrayObject *x_array, xparrayObject *y_array)
{
	unsigned i, j;

	PyObject *ret_val;
	/*Returns a reference to ret_val*/
	ret_val = PyObject_CallFunctionObjArgs(int_cb_callable, (PyObject *) x_array, (PyObject *) y_array, NULL);

	if (ret_val && !PyXParray_Check(ret_val))
	{
	   	xparrayObject *ret_array = (xparrayObject *) ret_val;
		/*vec = ret_array->data[0];*/ /*vec now points to ret_array->data[0]*/
		for (i = 0; i < 14; i++)
		{
			for (j = 0; j < 14; j++)
			{
				z_array->data[i][j] = ret_array->data[i][j];
			}
		}
	}
	else 
	{
		PyErr_SetString(PyExc_TypeError, "Callable object int2: bad behavior");
	}

	Py_XDECREF(ret_val);
}

/* Python interface of integral2 */

static PyObject *PYintegral2(PyObject *self, PyObject *args, PyObject *kwargs)
{
	/*integral arguments*/
	double AbsTol = 1e-5;
	double RelTol = 0;
	int Sector = 0;
	int Singular = 0;
	PyObject *limits = NULL;
	PyObject *py_Sector = NULL;
	PyObject *py_Singular = NULL;
	/*Values to be returned*/
	double error = 0;
	double result;
	PyObject *ret_tuple = PyTuple_New(2); 
	PyObject *ret_integral;
	PyObject *ret_error;

	static char *kwlist[] = {"f", "limits", "AbsTol", "RelTol", "Sector", "Singular", NULL};

	/* parse the input tuple with keywords */

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|$ddO!O!", kwlist, &int_cb_callable, &limits, &AbsTol, &RelTol, &PyBool_Type, &py_Sector, &PyBool_Type, &py_Singular))
	{
		return NULL;
	}

	if (py_Sector)
		Sector = PyObject_IsTrue(py_Sector);

	if (py_Singular)
		Singular = PyObject_IsTrue(py_Singular);

	/* ensure the first parameter is a callable */
	if (!PyCallable_Check(int_cb_callable))
	{
		PyErr_SetString(PyExc_TypeError, "integral: a callable is required\n");
		return NULL;
	}

	/* Checks if the limits are correct */
	if (!PyList_CheckExact(limits))
	{
		PyErr_SetString(PyExc_TypeError, "The argument limits only accept a list of lists\n");
		return NULL;
	}
	
	/* Ok, we have a list. Now checks if the list size is ok */
	if (PyList_Size(limits) != 2)
	{
		PyErr_SetString(PyExc_TypeError, "The argument limits only accept a list of size 2\n");
		return NULL;
	}

	/* The size of the list is ok. Now checks if the first element is a list of size 2 with
	 * two elements of type double
	 */
	PyObject *limits_ab = PyList_GetItem(limits, 0);
	if (!PyList_CheckExact(limits_ab))
	{
		PyErr_SetString(PyExc_TypeError, "The first element of limits must be a list\n");
		return NULL;
	}

	/* Ok, we have a list. Now checks if the list size is ok */
	if (PyList_Size(limits_ab) != 2)
	{
		PyErr_SetString(PyExc_TypeError, "The first element of limits must be a list of size 2\n");
		return NULL;
	}

	/*Now, checks if the two elements of limits_ab are of type double*/
	PyObject *py_a = PyList_GetItem(limits_ab, 0);
	PyObject *py_b = PyList_GetItem(limits_ab, 1);

	if (!PyNumber_Check(py_a) || !PyNumber_Check(py_b))
	{
		PyErr_SetString(PyExc_TypeError, "The first element of limits must be a list of size 2 with numbers\n");
		return NULL;
	}

	double a = PyFloat_AsDouble(py_a);
	double b = PyFloat_AsDouble(py_b);

	/* Now checks if the second element is a list of size 2 with
	 * two elements of type double or callable
	 */
	PyObject *limits_cd = PyList_GetItem(limits, 1);
	if (!PyList_CheckExact(limits_cd))
	{
		PyErr_SetString(PyExc_TypeError, "The second element of limits must be a list\n");
		return NULL;
	}

	/* Ok, we have a list. Now checks if the list size is ok */
	if (PyList_Size(limits_cd) != 2)
	{
		PyErr_SetString(PyExc_TypeError, "The second element of limits must be a list of size 2\n");
		return NULL;
	}

	/*Now, checks if the two elements of limits_cd are of type double or callable*/
	PyObject *py_c = PyList_GetItem(limits_cd, 0);
	PyObject *py_d = PyList_GetItem(limits_cd, 1);

	/*Reset global variables*/
	d_scalar = c_scalar = 0;

	if (PyNumber_Check(py_c))
	{
		c_scalar = PyFloat_AsDouble(py_c);
		c_callable = 0;
	} 
	else
	{
		if (!PyCallable_Check(py_c))
		{
			PyErr_SetString(PyExc_TypeError, "integral2: The lower limit of the first integral must be a number or a callable\n");
			return NULL;
		}
		c_callable = 1;
	}

	if (PyNumber_Check(py_d))
	{
		d_scalar = PyFloat_AsDouble(py_d);
		d_callable = 0;
	} 
	else
	{
		if (!PyCallable_Check(py_d))
		{
			PyErr_SetString(PyExc_TypeError, "integral2: The upper limit of the first integral must be a number or a callable\n");
			return NULL;
		}
		d_callable = 1;
	}

	/* Calls the integral function */
	result = integral2(int2_callback, a, b, py_c, py_d, AbsTol, RelTol, Sector, Singular, &error);

	ret_integral = Py_BuildValue("d", result);
	ret_error = Py_BuildValue("d", error);

	PyTuple_SetItem(ret_tuple, 0, ret_integral); /*steals the reference of ret_integral*/
	PyTuple_SetItem(ret_tuple, 1, ret_error); /*steals the reference of ret_error*/

	return ret_tuple;
}

/*****************************************************************************
 * 
 * 							Top level function
 * 
 *****************************************************************************
 */

PyObject *PyXP_Integral(PyObject *self, PyObject *args, PyObject *kwargs)
{
	int_cb_callable = PyTuple_GetItem(args, 0); /* Borrowed reference */

	if (!PyCallable_Check(int_cb_callable))
	{
		PyErr_SetString(PyExc_TypeError, "integral: a callable is required\n");
		return NULL;
	}

	PyObject *fc = PyObject_GetAttrString(int_cb_callable, "__code__");
	if (fc)
	{
		unsigned count = 0;
		PyObject *ac = PyObject_GetAttrString(fc, "co_argcount");
		if (ac)
		{
			count = PyLong_AsLong(ac);
			Py_DECREF(ac);
		}

		switch (count)
		{
			case 1:
				return PYintegral1(self, args, kwargs);
			case 2:
				return PYintegral2(self, args, kwargs);
			default:
				; /* nothing */		
		}

		Py_DECREF(fc);
	}

	return NULL;
}