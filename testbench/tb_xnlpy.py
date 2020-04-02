import xnlpy as xp
import math as m
from timeit import default_timer as timer
from scipy.integrate import quad

phi = (1 + m.sqrt(5)) / 2
acot = lambda x: m.pi/2 - m.atan(x)
true_res = 4 * m.pi * acot(m.sqrt(phi))

a = 0
b = m.pi

############# XNLPY #############

print("Numerical integration from xnlpy...")

func = lambda x: m.tan(x) / m.tan(x/2) * m.log((2*m.cos(x)**2 + 2*m.cos(x) + 1) / (2*m.cos(x)**2 - 2*m.cos(x) + 1))

init = timer()

result, error = xp.integral(func, a, b, points=15)

end = timer()

dt = (end - init)*10**6 # in us

print("Result:", result,"\nEst. Error:", error,"\nTrue Error:", result - true_res)

print("Elapsed time (in us): ",dt)

############# SCIPY #############

print("Numerical integration from scipy...")

init = timer()

result, error = quad(func, a, b)

end = timer()

dt = (end - init)*10**6 # in us

print("Result:", result,"\nEst. Error:", error,"\nTrue Error:", result - true_res)

print("Elapsed time (in us): ",dt)