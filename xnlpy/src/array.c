/***************************************************************************
XNL: eXtended Numerical Library
2020, Gilberto Jose Guimaraes de Sousa Mourao

XNL is a library made by some students of UFMA (Universidade Federal do 
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

static int ParseArguments(unsigned long *arr, Py_ssize_t size, PyObject *args)
{
	Py_ssize_t i;
	PyObject *temp_p, *temp_p2;

	for (i = 0; i < size; i++)
	{
		temp_p = PyTuple_GetItem(args, i);
		if (temp_p == NULL)
			return 0;

		/*Check if temp_p is numeric*/
		if (PyNumber_Check(temp_p) != 1)
		{
			PyErr_SetString(PyExc_TypeError, "Non-numeric argument.");
			return 0;
		}

		/*Convert number to python long and then C unsigned long*/
		temp_p2 = PyNumber_Long(temp_p);
		arr[i] = PyLong_AsUnsignedLong(temp_p2);
		Py_DECREF(temp_p2);
		if (arr[i] <= 0)
		{
			PyErr_SetString(PyExc_TypeError, "Non-positive number doesn't allowed as argument");
			return 0;
		}

		if (PyErr_Occurred())
			return 0;
	}

	return 1;
}

/*Creates an array with zeros*/
PyObject *py_zeros(PyObject *self, PyObject *args)
{
	Py_ssize_t TupleSize = PyTuple_Size(args);
	Py_ssize_t i;
	unsigned long *nums;
	PyObject *main_list;

	if (!TupleSize)
	{
		if (!PyErr_Occurred())
			PyErr_SetString(PyExc_TypeError, "You must supply at least one argument.");
		return NULL;
	}

	if (TupleSize > 2)
	{
		PyErr_SetString(PyExc_TypeError, "Currently, only 1d and 2d arrays are supported.");
		return NULL;
	}

	nums = malloc(TupleSize * sizeof(unsigned long));

	if(!ParseArguments(nums, TupleSize, args)){
        /* Make a cleanup and then return null*/
        return NULL;
    }

    /*Valid arguments. Continue...*/
    main_list = PyList_New(0);

    Py_ssize_t rows = TupleSize > 1 ? nums[0] : TupleSize;
    unsigned long cols = TupleSize > 1 ? nums[1] : nums[0];

    /*for loop to create the array*/
    for (i = 0; i < rows; i++)
    {
    	PyObject *list = PyList_New(cols);
    	unsigned long j;

    	for(j = 0; j < cols; j++)
    	{
    		PyObject *zero = PyLong_FromLong(0);
    		if(PyList_SetItem(list, j, zero))
    		{
    			PyErr_SetString(PyExc_TypeError, "Failed to create the array.");
    			return NULL;
    		}
    	}

    	if (PyList_Append(main_list, list))
    	{
    		PyErr_SetString(PyExc_TypeError, "Failed to create the array.");
    		return NULL;
    	}
    }

    /*free memory*/
    free(nums);

    return main_list;
}