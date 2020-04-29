import xnlpy as xp
from timeit import default_timer as timer

A = xp.array([1, 2, 3], [4, 5, 6])

print('~Testing print...\n')
A.print(2)

print('\n~Testing read and write...\n')
print('~Read method:',A.read(1,0))
A.write(1,0,99)
print('~After write:',A.read(1,0))
print('\n~Now printing...\n')
A.print(2)

print('\n~Testing zeros...\n')
B = xp.zeros(3,5)
B.print(2)

print('\n~Testing ones...\n')
C = xp.ones(5,3)
C.print(2)

print('\n~Testing eye...\n')
D = xp.eye(5)
D.print(2)
