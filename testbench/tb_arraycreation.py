import xnlpy as xp
#from timeit import default_timer as timer

print('\n~Testing zeros...\n')

B = xp.zeros(3,5)
print(B)

print('\n~Testing ones...\n')

C = xp.ones(5,3)
print(C)

print('\n~Testing eye...\n')

D = xp.eye(5)
print(D)

print('\n~Testing mult...\n')

mat1 = xp.array([1, 2, 3],[4, 5, 6]);
mat2 = xp.array([1, 2], [3, 4], [5, 6]);
mat1x2 = mat1 * mat2
print(mat1x2)
print('~1\n')
print(mat2 * 1.5)

print('\n~Testing transpose...\n')

xp.transpose(mat1).print(2)
test = xp.array([1, 2, 3, 4, 5])
print('~2\n')
print(test * xp.transpose(test))

print('\n~Testing plus and minus...\n')

matA = xp.array([1, 2], [3, 4], [5, 6])
matB = xp.array([1, 0], [0, 3], [2, 1])
print(matA + matB)
print('~3\n')
print(matA - matB)
print('~4\n')
print(matA + 2 + matB - 2)
