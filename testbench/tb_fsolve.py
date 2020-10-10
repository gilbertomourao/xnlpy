import xnlpy as xp 
import math as m

func = lambda x: m.sin(x)
dfunc = lambda x: m.cos(x)

sec_root = xp.fsolve(func, x = [m.pi-1,m.pi+1])
newton_root = xp.fsolve(func, df = dfunc, x = [m.pi-1])

print('Secant method:', sec_root,'\nNewton method:', newton_root)