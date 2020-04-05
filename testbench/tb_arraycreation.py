import xnlpy as xp
import time

start = time.time()

A = xp.zeros(5,2)

end = time.time()

dt = (end - start) * 1e6

print(A)

print('Elapsed time: ', dt)