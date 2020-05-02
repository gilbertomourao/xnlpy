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
	{"zeros", (PyCFunction) py_zeros, METH_VARARGS, "Zero array creation"},
	{"ones", (PyCFunction) py_ones, METH_VARARGS, "Unity array creation"},
	{"eye", (PyCFunction) py_eye, METH_VARARGS, "Identity matrix creation"},
	{"mult", (PyCFunction) py_mult, METH_VARARGS, "Array multiplication"},
	{"transpose", (PyCFunction) py_transpose, METH_VARARGS, "Array transpose"},
	{"add", (PyCFunction) py_add, METH_VARARGS, "Array addition"},
	{"sub", (PyCFunction) py_sub, METH_VARARGS, "Array subtraction"},
	{"integral", (PyCFunction) py_integral, METH_VARARGS | METH_KEYWORDS, "Numerical integration"},
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

	return m;
}