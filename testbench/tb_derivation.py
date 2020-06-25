import xnlpy as xp
import math as m 

func = lambda x: m.sin(x)
realvalue = m.cos(2.3)

print('Real value: ', realvalue)

print('\nTesting diff with default arguments:\n')

xnlvalue = xp.diff(func, 2.3)

print('XNL value: ',xnlvalue,
	  '\nTrue error: ',abs(realvalue - xnlvalue))

print('\nNow with a different step:\n')

xnlvalue = xp.diff(func, 2.3, step=1e-6)

print('XNL value: ',xnlvalue,
	  '\nTrue error: ',abs(realvalue - xnlvalue))

print('\nNow with a different step and more iterations:\n')

xnlvalue = xp.diff(func, 2.3, step=1e-6, iterations=10)

print('XNL value: ',xnlvalue,
	  '\nTrue error: ',abs(realvalue - xnlvalue))
