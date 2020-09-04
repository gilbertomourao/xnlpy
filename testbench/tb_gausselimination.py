import xnlpy as xp

back_A = xp.array([-3, 2, -1], [6, -6, 7], [3, -4, 4])
back_b = xp.array([-1], [-7], [-6])

up_A = xp.array([-3, 2, -1], [6, -6, 7], [3, -4, 4])
up_b = xp.array([-1], [-7], [-6])

jordan_A = xp.array([-3, 2, -1], [6, -6, 7], [3, -4, 4])
jordan_b = xp.array([-1], [-7], [-6])

# Backward subs
print('Backwards subs test:\n')
back_x = xp.GaussElimination(back_A, back_b, solve=True)

back_A.print(2)
back_b.print(2)
back_x.print(2)

# Upwards subs
print('\nUpwards subs test:\n')
up_x = xp.GaussElimination(up_A, up_b, solve=True, direction="tril")

up_A.print(2)
up_b.print(2)
up_x.print(2)

# Gauss-Jordan
print('\nGauss-Jordan test:\n')
xp.GaussElimination(jordan_A, jordan_b, direction="tril")
xp.GaussElimination(jordan_A, jordan_b)

jordan_A.print(2)
jordan_b.print(2)