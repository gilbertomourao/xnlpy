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

static PyMethodDef xnl_methods[] = 
{
	{"acos", (PyCFunction) PyXP_Acos, METH_VARARGS, "acos computation"},
	{"acosh", (PyCFunction) PyXP_Acosh, METH_VARARGS, "acosh computation"},
	{"asin", (PyCFunction) PyXP_Asin, METH_VARARGS, "asin computation"},
	{"asinh", (PyCFunction) PyXP_Asinh, METH_VARARGS, "asinh computation"},
	{"atan", (PyCFunction) PyXP_Atan, METH_VARARGS, "atan computation"},
	{"atanh", (PyCFunction) PyXP_Atanh, METH_VARARGS, "atanh computation"},
	{"cos", (PyCFunction) PyXP_Cos, METH_VARARGS, "cos computation"},
	{"cosh", (PyCFunction) PyXP_Cosh, METH_VARARGS, "cosh computation"},
	{"sin", (PyCFunction) PyXP_Sin, METH_VARARGS, "sin computation"},
	{"sinh", (PyCFunction) PyXP_Sinh, METH_VARARGS, "sinh computation"},
	{"tan", (PyCFunction) PyXP_Tan, METH_VARARGS, "tan computation"},
	{"tanh", (PyCFunction) PyXP_Tanh, METH_VARARGS, "tanh computation"},
	{"exp", (PyCFunction) PyXP_Exp, METH_VARARGS, "exp computation"},
	{"log", (PyCFunction) PyXP_Log, METH_VARARGS, "log computation"},
	{"log10", (PyCFunction) PyXP_Log10, METH_VARARGS, "log10 computation"},
	{"sqrt", (PyCFunction) PyXP_Sqrt, METH_VARARGS, "sqrt computation"},
	{"ceil", (PyCFunction) PyXP_Ceil, METH_VARARGS, "ceil computation"},
	{"fabs", (PyCFunction) PyXP_Fabs, METH_VARARGS, "fabs computation"},
	{"floor", (PyCFunction) PyXP_Floor, METH_VARARGS, "floor computation"},
	

	{"zeros", (PyCFunction) PyXP_Zeros, METH_VARARGS, "Zero array creation"},
	{"ones", (PyCFunction) PyXP_Ones, METH_VARARGS, "Unity array creation"},
	{"eye", (PyCFunction) PyXP_Eye, METH_VARARGS, "Identity matrix creation"},
	{"transpose", (PyCFunction) PyXP_Transpose, METH_VARARGS, "Array transpose"},
	{"GaussElimination", (PyCFunction) PyXP_GaussElimination, METH_VARARGS | METH_KEYWORDS, "Gauss Elimination"},
	{"arange", (PyCFunction) PyXP_Arange, METH_VARARGS, "Range array creation"},

	{"fsolve", (PyCFunction) PyXP_Fsolve, METH_VARARGS | METH_KEYWORDS, "Numerical root finder for nonlinear function"},
	
	
	{"diff", (PyCFunction) PyXP_Diff, METH_VARARGS | METH_KEYWORDS, "Numerical differentiation"},
	{"integral", (PyCFunction) PyXP_Integral, METH_VARARGS | METH_KEYWORDS, "Numerical integration"},

	{"lambertw", (PyCFunction) PyXP_LambertW, METH_VARARGS, "Principal branch of Lambert W function"},
	
	{NULL, NULL, 0, NULL} /*Sentinel*/
};

static struct PyModuleDef xnl_module = 
{
	PyModuleDef_HEAD_INIT,
	"xnlpy", /*name of the module*/
	NULL, /*module documentation*/
	-1, /*still can't understand this*/

	xnl_methods
};

PyMODINIT_FUNC PyInit_xnlpy(void)
{
	PyObject *m;

	if (PyType_Ready(&xparrayType) < 0)
		return NULL;

	m = PyModule_Create(&xnl_module);

	if (m == NULL)
		return NULL;

	Py_INCREF(&xparrayType);
	if (PyModule_AddObject(m, "array", (PyObject *) &xparrayType) < 0)
	{
		Py_DECREF(&xparrayType);
		Py_DECREF(m);
		return NULL;
	}

	PyObject *pi = Py_BuildValue("d", XNL_PI);

	if (PyModule_AddObject(m, "pi", pi))
	{
		Py_DECREF(pi);
		Py_DECREF(m);
		return NULL;
	}

	PyObject *eps = Py_BuildValue("d", XNL_EPS);

	if (PyModule_AddObject(m, "eps", eps))
	{
		Py_DECREF(eps);
		Py_DECREF(m);
		return NULL;
	}

	PyObject *inf = Py_BuildValue("d", XNL_INFINITY);

	if (PyModule_AddObject(m, "inf", inf))
	{
		Py_DECREF(inf);
		Py_DECREF(m);
		return NULL;
	}

	PyObject *nan = Py_BuildValue("d", XNL_NAN);

	if (PyModule_AddObject(m, "nan", nan))
	{
		Py_DECREF(nan);
		Py_DECREF(m);
		return NULL;
	}

	return m;
}