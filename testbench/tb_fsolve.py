import xnlpy as xp

func = lambda x: xp.sin(x)
dfunc = lambda x: xp.cos(x)

sec_root = xp.fsolve(func, x = [xp.pi-1,xp.pi+1])
newton_root = xp.fsolve(func, df = dfunc, x = [xp.pi-1])

print('Secant method:', sec_root,'\nNewton method:', newton_root)