import xnlpy as xp
from timeit import default_timer as timer

start = timer()

A = xp.zeros(5,2)

end = timer()

dt = (end - start) * 1e6

print(A)

print('Elapsed time (in us): ', dt)