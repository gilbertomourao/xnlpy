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

#ifndef XNLPY_H
#define XNLPY_H

#include <sys/time.h> /*struct timeval declaration*/
#define PY_SSIZE_T_CLEAN
#include <Python.h>

/*Constants*/
#define XNL_PI 3.14159265358979323846
#define XNL_EPS 2.220446049250313e-16
#define XNL_INFINITY HUGE_VAL
#define XNL_NAN (0.0f / 0.0f)

/*Basic Math*/
PyObject *PyXP_Acos(PyObject *, PyObject *);
PyObject *PyXP_Acosh(PyObject *, PyObject *);
PyObject *PyXP_Asin(PyObject *, PyObject *);
PyObject *PyXP_Asinh(PyObject *, PyObject *);
PyObject *PyXP_Atan(PyObject *, PyObject *);
PyObject *PyXP_Atanh(PyObject *, PyObject *);
PyObject *PyXP_Cos(PyObject *, PyObject *);
PyObject *PyXP_Cosh(PyObject *, PyObject *);
PyObject *PyXP_Sin(PyObject *, PyObject *);
PyObject *PyXP_Sinh(PyObject *, PyObject *);
PyObject *PyXP_Tan(PyObject *, PyObject *);
PyObject *PyXP_Tanh(PyObject *, PyObject *);
PyObject *PyXP_Exp(PyObject *, PyObject *);
PyObject *PyXP_Log(PyObject *, PyObject *);
PyObject *PyXP_Log10(PyObject *, PyObject *);
PyObject *PyXP_Sqrt(PyObject *, PyObject *);
PyObject *PyXP_Ceil(PyObject *, PyObject *);
PyObject *PyXP_Fabs(PyObject *, PyObject *);
PyObject *PyXP_Floor(PyObject *, PyObject *);


/*Linear Algebra*/
#include "arrayobject.h"

xparrayObject *PyXP_Zeros(PyObject *, PyObject *);
xparrayObject *PyXP_Ones(PyObject *, PyObject *);
xparrayObject *PyXP_Eye(PyObject *, PyObject *);
xparrayObject *PyXP_Transpose(PyObject *, PyObject *);
xparrayObject *PyXP_GaussElimination(PyObject *, PyObject *, PyObject *);

/*Nonlinear equations*/
PyObject *PyXP_Fsolve(PyObject *, PyObject *, PyObject *);

/*Single Variable Calculus*/
PyObject *PyXP_Diff(PyObject *, PyObject *, PyObject *);
PyObject *PyXP_Integral(PyObject *, PyObject *, PyObject *);

#endif /*XNLPY_H*/