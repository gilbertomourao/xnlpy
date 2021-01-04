"""
Testing TwoD algorithm
"""
from xnlpy import *

problem = 1
fevals = 0

# First test
def f1(x, y):

	global fevals
	fevals = fevals + 1
	
	if (problem == 1):
		return (y**2 * x**2 + y**2 + x**2) * cos(x)

	if (problem == 2):
		return y**2 * sin(y + x)**2 * cos(x)

	if (problem == 3):
		return exp(sin(y)**2 * sin(x)**2) * cos(x)

	return 0

print("SSex\n")

problem = 1
fevals = 0

true_result = pi**5 / 3 - pi**3 / 3 - 8*pi

result, error = integral(f1, [[-pi/2, pi/2], [-pi, pi]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

problem = 2
fevals = 0

true_result = (2/3) * pi**3 - pi / 3

result, error = integral(f1, [[-pi/2, pi/2], [-pi, pi]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

problem = 3
fevals = 0

true_result = 15.24409680

result, error = integral(f1, [[-pi/2, pi/2], [-pi, pi]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# Second Test

def f2(x, y):

	global fevals
	fevals = fevals + 1

	if problem == 1:
		return exp(x * y)

	if problem == 2:
		return sin(pi/2 * (x + 2*y))

	if problem == 3:
		return 1 / ((x+1)**2 + (y+2)**2)

	return 0

print("RRex\n")

problem = 1
fevals = 0

true_result = 4.229003501502914

result, error = integral(f2, [[-1, 1], [-1, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

problem = 2
fevals = 0

true_result = 0

result, error = integral(f2, [[-1, 1], [-1, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

problem = 3
fevals = 0

true_result = 0.9379490574117436

result, error = integral(f2, [[-1, 1], [-1, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# Third Test

def f3(x, y):

	global fevals
	fevals = fevals + 1

	return 2*x*cos(y)

print("NRex\n")

fevals = 0

true_result = (-cos(9) - 9/2) - (-cos(1) - 1/2)

result, error = integral(f3, [[1, 3], [pi/6, lambda x: x*x]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# Fourth Test

def f4(x, y):

	global fevals
	fevals = fevals + 1

	return x + y

print("NAGex\n")

fevals = 0

true_result = 2/3

result, error = integral(f4, [[0, 1], [0, lambda x: sqrt(1 - x*x)]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# Singular set to 1

fevals = 0

result, error = integral(f4, [[0, 1], [0, lambda x: sqrt(1 - x*x)]], Singular = True)

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# Describe the region as a generalized sector

fevals = 0

result, error = integral(f4, [[0, pi/2], [0, 1]], Sector = True)

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# Fifth Test

def f5(x, y):

	global fevals
	fevals = fevals + 1

	if problem == 1:
		return cos(x + y)
	if problem == 2:
		return exp(fabs(x+y-1))
	if problem == 3:
		return exp(sin(x) * cos(y))
	if problem == 4:
		rsq = x*x + y*y
		return (x * rsq**1.5) / (rsq + 1e-6) - 100*rsq

	return 0

print("KahanerEx\n")

problem = 1
fevals = 0

true_result = 2*cos(1)-cos(2)-1

result, error = integral(f5, [[0, 1], [0, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

problem = 2
fevals = 0

true_result = 2*exp(1) - 4

result, error = integral(f5, [[0, 1], [0, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

problem = 3
fevals = 0

true_result = 1.508588079098596

result, error = integral(f5, [[0, 1], [0, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

problem = 4
fevals = 0

true_result = -800/3

result, error = integral(f5, [[-1, 1], [-1, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# Sixth Test

def f6(x, y):

	global fevals
	fevals = fevals + 1

	if problem == 1:
		return fabs(x*x + y*y - 1/4)
	if problem == 2:
		return 1.0 / (1 - x*y)
	if problem == 3:
		return sqrt(fabs(x - y))

	return 0

print("DRex\n")

problem = 1
fevals = 0

true_result = 5/3 + pi/16

result, error = integral(f6, [[-1, 1], [-1, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

problem = 2
fevals = 0

true_result = pi * pi / 6

result, error = integral(f6, [[0, 1], [0, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# With Singular set to true

fevals = 0
result, error = integral(f6, [[0, 1], [0, 1]], Singular = True)

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

problem = 3
fevals = 0

true_result = 8 / 15

result, error = integral(f6, [[0, 1], [0, 1]])

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# Twice the integral over [0,1]x[0,x]

fevals = 0
result, error = integral(f6, [[0, 1], [0, lambda x: x]])
result *= 2
error *= 2

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# With Sungular set to True

fevals = 0
result, error = integral(f6, [[0, 1], [0, lambda x: x]], Singular = True)
result *= 2
error *= 2

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

# Seventh Test

def ft(x, y):

	global fevals
	fevals = fevals + 1

	return exp(-1*(x*x + y*y))

print("Gaussian\n")

fevals = 0

true_result = pi
result, error = integral(ft, [[0,inf], [0,inf]], Singular = True)
result *= 4
error *= 4

print('Result = ',result,'\nError = ',error,'\nTrue error = ', true_result - result,'\nFevals = ', fevals,'\n')

