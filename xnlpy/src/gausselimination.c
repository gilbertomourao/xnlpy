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
#include <string.h>
#include <inttypes.h>
#include <stdbool.h>

#define FLAG(I, N, DIR) ((DIR == 1) ? (I < N) : (I >= 0))

typedef struct Complete_Array
{
	double **data;
	unsigned rows;
	unsigned cols;
} Array;

/**
 * Functions
 */

static void destroyArray(Array *my_arr)
{
	int i;

	if (!(my_arr && my_arr->rows && my_arr->cols))
	{
		printf("ERROR: In destroyArray. Array does not have allocated memory.\n");
		return;
	}

	for (i = 0; i < my_arr->rows; i++)
	{
		free(my_arr->data[i]);
		my_arr->data[i] = NULL;
	}

	free(my_arr->data);
	my_arr->data = NULL;

	free(my_arr);
	my_arr = NULL;
}

static void normalize(Array *mat_A, Array *vec_b)
{
	int i, j;
	double data, absdata, max_val;

	for (i = 0; i < mat_A->rows; i++)
	{
		max_val = 0;

		for (j = 0; j < mat_A->cols; j++)
		{
			data = mat_A->data[i][j];
			absdata = (data >= 0) ? data : (-data);
			if (absdata > max_val)
			{
				max_val = absdata;
			}
		}

		if (max_val > 0)
		{
			for (j = 0; j < mat_A->cols; j++)
			{
				mat_A->data[i][j] /= max_val;
			}
			vec_b->data[i][0] /= max_val;
		}
	}
}

static void switch_rows(Array *mat_A, Array *vec_b, unsigned position, unsigned mrow)
{
	int j;
	double aux;

	for (j = position; j < mat_A->cols; j++)
	{
		aux = mat_A->data[position][j];
		mat_A->data[position][j] = mat_A->data[mrow][j];
		mat_A->data[mrow][j] = aux;
	}

	aux = vec_b->data[position][0];
	vec_b->data[position][0] = vec_b->data[mrow][0];
	vec_b->data[mrow][0] = aux;
}

static void switch_cols(Array *mat_A, unsigned position, unsigned mcol)
{
	int i;
	double aux;

	for (i = position; i < mat_A->rows; i++)
	{
		aux = mat_A->data[i][position];
		mat_A->data[i][position] = mat_A->data[i][mcol];
		mat_A->data[i][mcol] = aux;
	}
}

static int row_pivoting(Array *mat_A, Array *vec_b, unsigned position)
{
	int i;
	double data, absdata, max_val = 0; 
	unsigned mrow = position;

	/*Checks for the max value at current column*/
	for (i = position; i < mat_A->rows; i++)
	{
		data = mat_A->data[i][position];
		absdata = (data >= 0) ? data : (-data);
		if (absdata > max_val)
		{
			max_val = absdata;
			mrow = i;
		}
	}

	/*In success, performs switch operation and return 1. If failed, return 0.*/
	if (max_val != 0)
	{
		switch_rows(mat_A, vec_b, position, mrow);
		return 1;
	} else
	{
		return 0;
	}
}

static int col_pivoting(Array *mat_A, unsigned position)
{
	int j;
	double data, absdata, max_val = 0; 
	unsigned mcol = position;

	/*Checks for the max value at current row*/
	for (j = position; j < mat_A->cols; j++)
	{
		data = mat_A->data[position][j];
		absdata = (data >= 0) ? data : (-data);
		if (absdata > max_val)
		{
			max_val = absdata;
			mcol = j;
		}
	}

	/*In success, performs switch operation and return 1. If failed, return 0.*/
	if (max_val != 0)
	{
		switch_cols(mat_A, position, mcol);
		return 1;
	} else
	{
		return 0;
	}
}

static int complete_pivoting(Array *mat_A, Array *vec_b, unsigned position)
{
	int i, j;
	double data, absdata, max_val = 0; 
	unsigned mrow = position, mcol = position;

	/*Checks for the max value in the sub-matrix*/
	for (i = position; i < mat_A->rows; i++)
	{
		for (j = position; j < mat_A->cols; j++)
		{
			data = mat_A->data[i][j];
			absdata = (data >= 0) ? data : (-data);
			if (absdata > max_val)
			{
				max_val = absdata;
				mrow = i;
				mcol = j;
			}
		}
	}

	/*In success, performs switch operation and return 1. If failed, return 0.*/
	if (max_val != 0)
	{
		switch_rows(mat_A, vec_b, position, mrow);
		switch_cols(mat_A, position, mcol);
		return 1;
	} else
	{
		return 0;
	}
}

static double ge_summation(Array *vec_sol, Array *mat_A, int start, int end)
{
	int k;
	double result = 0;
	unsigned position = (start == 0) ? (end + 1) : (start - 1); 

	for (k = start; k <= end;k++)
	{
		result += vec_sol->data[k][0] * mat_A->data[position][k];
	}

	return result;
}

/*xnl function*/
static Array *xnl_GaussElimination(Array *mat_A, Array *vec_b, char *direction, char *pivot, bool normal, bool solve)
{
	int h, i, j;
	double factor, norm_val;
	unsigned n;
	uint8_t pivoting_type;
	int8_t dir;

	/*chek if mat_A and vec_b are valid parameters*/
	if (!mat_A)
	{
		printf("WARNING: In GaussElimination. First argument cannot be a NULL pointer and must "
			   "be specified. No operation will be executed.\n");
		return NULL;
	}

	/*check if mat_A is a square matrix*/
	if (mat_A->rows != mat_A->cols)
	{
		printf("WARNING: In GaussElimination. First argument must be a square matrix. " 
				"No operation will be executed.\n");
		return NULL;
	}

	n = mat_A->rows;

	if (!vec_b)
	{
		printf("WARNING: In GaussElimination. Second argument cannot be a NULL pointer and must "
			   "be specified. No operation will be executed.\n");
		return NULL;
	}

	/*check if b is a 1D vector and length(b) = length(A)*/
	if (vec_b->rows != n || vec_b->cols != 1)
	{
		printf("WARNING: In GaussElimination. The second argument must be a 1D vector " 
				"with the same length of the first argument. No operation will be executed.\n");
		return NULL;
	}

	/*check if direction is a valid string*/
	if (!strcmp(direction, "tril")) dir = -1;
	else if (!strcmp(direction, "triu")) dir = 1;
	else
	{
		printf("ERROR: In GaussElimination. Invalid parameter direction. The valid parameters are tril and triu. "
			   "No operation will be executed.\n");
		return NULL;
	}

	/*check if pivot is a valid string*/
	if (!pivot) pivoting_type = 0;
	else if (!strcmp(pivot, "row")) pivoting_type = 1;
	else if (!strcmp(pivot, "col")) pivoting_type = 2;
	else if (!strcmp(pivot, "complete")) pivoting_type = 3;
	else 
	{
		printf("ERROR: In GaussElimination. Invalid parameter pivot. The valid parameters are row, col and complete. " 
			   "No operation will be executed.\n");
		return NULL;
	}

	/*check if the user wants to normalize [A|b]*/
	if (normal) normalize(mat_A, vec_b);

	for (h = (n-1)*((1-dir)>>1); FLAG(h, n, dir); h+=dir)
	{
		switch (pivoting_type)
		{
			case 1: 
			{
				if (!row_pivoting(mat_A, vec_b, h))
				{
					printf("WARNING: In GaussElimination. Pivoting method failed. Only zeros in row %d.\n", h);
					return NULL;
				}
				break;
			}
			case 2:
			{
				if (!col_pivoting(mat_A, h))
				{
					printf("WARNING: In GaussElimination. Pivoting method failed. Only zeros in column %d.\n", h);
					return NULL;
				}
				break;
			}
			case 3:
			{
				if (!complete_pivoting(mat_A, vec_b, h))
				{
					printf("WARNING: In GaussElimination. Pivoting method failed. Only zeros in submatrix.\n");
					return NULL;
				}
				break;
			}
			default: ;
		} /*end switch pivoting_type*/

		/*Normalize relative to the main diagonal*/
		norm_val = mat_A->data[h][h];
		for (j = 0; j < n; j++)
		{
			mat_A->data[h][j] /= norm_val;
		}
		vec_b->data[h][0] /= norm_val;

		for (i = h+dir; FLAG(i, n, dir); i+=dir)
		{
			factor = mat_A->data[i][h];
			for (j = h*((1+dir)>>1); FLAG(h-j,(int)(n+h-(j<<1)),dir); j++)
			{
				mat_A->data[i][j] -= factor * mat_A->data[h][j];
			}
			vec_b->data[i][0] -= factor * vec_b->data[h][0];
		}
	}

/*
	printArray(mat_A, 2);
	printArray(vec_b, 2);
*/

	if (solve == true) ; /*do nothing, just continue to the solution*/
	else if (solve != false)
	{
		printf("ERROR: In GaussElimination. Invalid parameter solve. The valid parameters are solve and NULL. " 
			   "No operation will be executed.\n");
		return NULL;
	} else return NULL;

	/*finally, performs top/backward substitution to find the solution*/

	/*First, create the solution array*/
	Array *vec_sol;

	vec_sol = malloc(sizeof(Array));

	if (vec_sol == NULL)
	{
		printf("ERROR: In GaussElimination. Memory allocation failed! During the creation of solution vector.\n");
		exit(EXIT_FAILURE);
	}

	vec_sol->data = malloc(n * sizeof(double *));

	if (vec_sol->data == NULL)
	{
		printf("ERROR: In GaussElimination. Memory allocation failed! During the creation of solution vector's data.\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < n; i++)
	{
		vec_sol->data[i] = malloc(sizeof(double));
		if (vec_sol->data[i] == NULL)
		{
			printf("ERROR: In GaussElimination. Memory allocation failed! During the creation of solution vector's data[%d].\n", i);
			exit(EXIT_FAILURE);
		}
	}

	vec_sol->rows = n;
	vec_sol->cols = 1;

	for (i = (n-1)*((dir + 1) >> 1); FLAG(i, n, -dir); i-=dir)
	{
		vec_sol->data[i][0] = vec_b->data[i][0] - ge_summation(vec_sol, mat_A, (i+1)*((1+dir)>>1), ((n+i + (n-i)*dir) >> 1) - 1);
	}

	return vec_sol;
}

/**
 * Python interface
 */

xparrayObject *PyXP_GaussElimination(PyObject *self, PyObject *args, PyObject *kwargs)
{
	/*arguments*/
	PyObject *arg_mat_A = NULL, *arg_vec_b = NULL;
	char *direction = "triu", *pivot = NULL;
	PyObject *py_normal = NULL, *py_solve = NULL; 
	bool normal = false, solve = false;

	/*keyword names*/
	static char *kwlist[] = {"A", "b", "direction", "pivot", "normal", "solve", NULL};

	/*input validation*/
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|$ssO!O!", kwlist, &arg_mat_A, &arg_vec_b, &direction, &pivot, &PyBool_Type, &py_normal, &PyBool_Type, &py_solve))
	{
		return NULL;
	}

	if (py_normal)
	{
		normal = PyObject_IsTrue(py_normal);
	}

	if (py_solve)
	{
		solve = PyObject_IsTrue(py_solve);
	}

	/*xnlpy.array arguments*/
	/***************************************************************************************/
	if (PyXParray_Check(arg_mat_A))
	{
		PyErr_SetString(PyExc_TypeError, "The first argument's type must be xnlpy.array.");
		return NULL;
	}

	if (PyXParray_Check(arg_vec_b))
	{
		PyErr_SetString(PyExc_TypeError, "The second argument's type must be xnlpy.array.");
		return NULL;
	}

	/*Passed*/
	xparrayObject *array_A = (xparrayObject *)arg_mat_A,
				  *array_b = (xparrayObject *)arg_vec_b;

	Array mat_A, vec_b;

	mat_A.data = array_A->data;
	mat_A.rows = array_A->rows;
	mat_A.cols = array_A->cols;

	vec_b.data = array_b->data;
	vec_b.rows = array_b->rows;
	vec_b.cols = array_b->cols;	

	/*return object*/
	Array *vec_sol = xnl_GaussElimination(&mat_A, &vec_b, direction, pivot, (bool)normal, (bool)solve);

	if (vec_sol)
	{
		/*must construct the solution vector as a xnlpy.array object*/

		/*Initialize a new array*/
		PyObject *argList = PyTuple_New(vec_sol->rows);

		Py_ssize_t i, j;
		for (i = 0; i < vec_sol->rows; i++)
		{
			PyObject *List = PyList_New(vec_sol->cols);

			for (j = 0; j < vec_sol->cols; j++)
			{
				PyObject *temp = Py_BuildValue("d", vec_sol->data[i][j]);

				if (PyList_SetItem(List, j, temp) == -1) /*steals the reference of temp*/
				{
					return NULL;
				}
			}

			PyTuple_SetItem(argList, i, List); /*steals the reference of List*/
		}

		/*call the array object*/
		xparrayObject *array_sol = (xparrayObject *)PyObject_CallObject((PyObject *) &xparrayType, argList);

		Py_DECREF(argList);

		if (array_sol == NULL)
		{
			/*Will print the messages from xparray*/
			return NULL;
		}

		/*Destroys vec_sol*/
		destroyArray(vec_sol);

		/*Success on the creation, then return the solution array*/
		return array_sol;
	}
	else 
	{
		return (xparrayObject *)Py_BuildValue("");
	}
}