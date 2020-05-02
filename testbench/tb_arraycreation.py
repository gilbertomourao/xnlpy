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

print('\n~Testing mult...\n')
mat1 = xp.array([1, 2, 3],[4, 5, 6]);
mat2 = xp.array([1, 2], [3, 4], [5, 6]);
mat1x2 = xp.mult(mat1, mat2)
mat1x2.print(2)

print('\n~Testing transpose...\n')
xp.transpose(mat1).print(2)
test = xp.array([1, 2, 3, 4, 5])
print('\n')
xp.mult(test, xp.transpose(test)).print(2)

print('\n~Testing add and sub...\n')
matA = xp.array([1, 2], [3, 4], [5, 6])
matB = xp.array([1, 0], [0, 3], [2, 1])
xp.add(matA, matB).print(2)
xp.sub(matA, matB).print(2)