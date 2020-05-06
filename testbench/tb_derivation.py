import xnlpy as xp
import math as m 

func = lambda x: m.sin(x)

xnlvalue = xp.diff(func, 2.3)
realvalue = m.cos(2.3)

print('XNL value: ',xnlvalue,
	  '\nReal value: ',realvalue)