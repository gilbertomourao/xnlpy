import xnlpy as xp
import math as m
from timeit import default_timer as timer

############# Function 1 #############

print("Integrating the first function...")

phi = (1 + m.sqrt(5)) / 2
acot = lambda x: m.pi/2 - m.atan(x)
true_res = 4 * m.pi * acot(m.sqrt(phi))

a = 0
b = m.pi

func = lambda x: m.tan(x) / m.tan(x/2) * m.log((2*m.cos(x)**2 + 2*m.cos(x) + 1) / (2*m.cos(x)**2 - 2*m.cos(x) + 1))

init = timer()
result, error = xp.integral(func, a, b, points=15)
end = timer()

dt = (end - init)*10**6 # in us

print("Result:", result,"\nEst. Error:", error,"\nTrue Error:", result - true_res)
print("Elapsed time (in us): ",dt)

############# Function 2 #############

print("Integrating the second function...")

true_res = (m.pi / 8) * m.log(m.pi**2 / 8)

a = 0
b = 1

func = lambda x: m.atan( (m.atanh(x) - m.atan(x)) / (m.pi + m.atanh(x) - m.atan(x)) ) / x

init = timer()
result, error = xp.integral(func, a, b, points=10, tolerance=1e-15, depth=50)
end = timer()

dt = (end - init)*10**6 # in us

print("Result:", result,"\nEst. Error:", error,"\nTrue Error:", result - true_res)
print("Elapsed time (in us): ",dt)

############# Function 3 #############

print("Integrating the third function...")

true_res = m.pi / 2

a = 0
b = m.inf

func = lambda x: (x**8 - 4*x**6 + 9*x**4 - 5*x**2 + 1) / (x**12 - 10*x**10 + 37*x**8 - 42*x**6 + 26*x**4 - 8*x**2 + 1)

init = timer()
result, error = xp.integral(func, a, b, points=15, tolerance=1e-15)
end = timer()

dt = (end - init)*10**6 # in us

print("Result:", result,"\nEst. Error:", error,"\nTrue Error:", result - true_res)
print("Elapsed time (in us): ",dt)