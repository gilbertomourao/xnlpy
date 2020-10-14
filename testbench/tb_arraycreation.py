import xnlpy as xp
#from timeit import default_timer as timer

A = xp.array([1, 2, 3], [4, 5, 6])

print('~Testing print...\n')
A.print(2)

print('\n~Testing read and write...\n')
print('~Read method:',A[1][0])
print(len(A)) # returns the largest array dimension
A[1][0] = 99
print('\n~Now printing...\n')
A[0:1][0].print(2)

print('\n~test...\n')

test = xp.array([1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5])
test[4:3:-2].print(1)

print('\n~Testing xparray to xparray assignment...\n')

test[2] = xp.zeros(1,5)
test.print(1)

print(test[2][2:2])

print('\n~Testing xparray to xparray assignment...\n')

test[2:3] = xp.array([1,0,1,0,1],[1,0,1,0,1])
test.print(1)

print('\n~Testing number to xparray assignment...\n')

test[4:0:-3][1:4:2] = 99
test.print(1)

print('\n~Testing xparray to xparray assignment...\n')

test[4:0:-3][1:4:2] = xp.array([1,2],[1,2])
test.print(1)

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
mat1x2 = mat1 * mat2
mat1x2.print(2)
print('~1\n')
(mat2 * 1.5).print(2)

print('\n~Testing transpose...\n')

xp.transpose(mat1).print(2)
test = xp.array([1, 2, 3, 4, 5])
print('~2\n')
print(test * xp.transpose(test))

print('\n~Testing plus and minus...\n')

matA = xp.array([1, 2], [3, 4], [5, 6])
matB = xp.array([1, 0], [0, 3], [2, 1])
(matA + matB).print(2)
print('~3\n')
(matA - matB).print(2)
print('~4\n')
(matA + 2 + matB - 2).print(2)
