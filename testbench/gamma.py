import xnlpy as xp

pgamma = lambda x: xp.integral(lambda t, x: t**(x-1)*xp.exp(-t), 0, xp.inf, args=(x,), depth=50)

