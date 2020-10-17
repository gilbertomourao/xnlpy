import xnlpy as xp 

A = xp.array([1, 2, 3], [4, 5, 6])

print('~Testing print...\n')
print(A)

print('\n~Testing read and write...\n')
print('~Read method:',A[1][0])
print('~Length of A = %d'%len(A)) # returns the largest array dimension
A[1][0] = 99
print('\n~Now printing...\n')
print(A[0:1][0])

print('\n~test...\n')

test = xp.array([1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5])
print(test[4:3:-2])

print('\n~Testing xparray to xparray assignment...\n')

test[2] = xp.zeros(1,5)
print(test)

print(test[2][2:2])

print('\n~Testing xparray to xparray assignment...\n')

test[2:3] = xp.array([1,0,1,0,1],[1,0,1,0,1])
print(test)

print('\n~Testing number to xparray assignment...\n')

test[4:0:-3][1:4:2] = 99
print(test)

print('\n~Testing xparray to xparray assignment...\n')

test[4:0:-3][1:4:2] = xp.array([1,2],[-1,-2])
print(test)

"""

1D array indexing and assignment

"""

print('\n~1D array...\n')

A1d = xp.array([1,2,3,4,5,6,7,8,9])

print('\n~Testing read and write...\n')
print('~Read method:',A1d[1])
print('~Length of A = %d'%len(A1d)) # returns the largest array dimension
A1d[1] = 99
print('\n~Now printing...\n')
print(A1d[0:1])