import xnlpy as xp
import time

def lambertw(x):
	"""
	Calcula o valor de W(x) para x real. Não 
	trabalha no domínio complexo. 

	A representação de W(x) como uma integral foi dada por Istvàn Mezö.

	Ver: https://arxiv.org/pdf/2012.02480.pdf

	O valor máximo suportado pela função fica entre 0.661e308 e 0.662e308.
	Os valores dos parâmetros tolerance e depth foram obtidos empiricamente.
	"""
	f_int = lambda t, k: xp.log(1 + k*xp.sin(t)/t*xp.exp(t/xp.tan(t)))

	return 1/xp.pi * xp.integral(f_int, 0, xp.pi, args = (x,), tolerance = 1e-15, depth = 30)[0]

if __name__ == '__main__':

	# Limite mínimo entre -0.366 e -0.365
	print('W(-0.365) =',lambertw(-0.365))

	# Limite máximo entre 0.661e308 e 0.662e308
	print('W(0.661e308) =',lambertw(0.661e308))

	# Teste de tempo para calcular a perda de carga h

	L = 5 # comp em m
	Q = 2.35 * 1e-3 # L/s -> m³/s
	DI = 20 * 1e-3 # m
	A = xp.pi*DI**2/4 # m²
	V = Q/A # m/s

	P = xp.pi*DI # perímetro molhado em m
	Rh = A / P # Raio hidráulico em m

	nu = 1.003e-6 # m²/s 

	Re = V*DI/nu # Reynolds
	rug = 0.015e-3 # rugosidade em m

	a = 2.51 / Re
	b = rug / (14.8 * Rh)

	start = time.time()
	ln_10 = xp.log(10)
	f = (2*lambertw(ln_10/(2*a)*10**(b/(2*a)))/ln_10 - b/a)**-2

	print(f)

	g = 9.8 # m²/s

	h = f*(L/DI)*V**2/(2*g)
	dt = time.time() - start

	print(h, dt)

	# Utilizando aproximação de Fair-Whipple-Hsiao

	start = time.time()

	h = 8.69e-4*(Q**1.75)*(DI**-4.75)*L
	dt = time.time() - start

	print(h, dt)