import xnlpy as xp 

f = lambda x, a, b: a*x + b

print(xp.integral(f, 0, 1, args=(1,1)))
print(xp.integral(f, 0, 1, args=(1,1,)))

#print(xp.integral(f,0,1, args=1))
#print(xp.integral(f,0,1, args=(1,)))
#print(xp.integral(f,0,1, args=[1,1]))
#print(xp.integral(f,0,1, args=(1,1,1)))
#print(xp.integral(f,0,1, args=(1,'a')))
#print(xp.integral(f,0,1))

a = 0
b = xp.pi

func = lambda x: x

print(xp.integral(func,0,1))