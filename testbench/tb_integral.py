import xnlpy as xp
#from timeit import default_timer as timer
#import statistics
import matplotlib.pyplot as plt

'''
def time_integral(iterations, func, a, b, **kwargs):
    points = kwargs.get('points')
    tolerance = kwargs.get('tolerance')
    depth = kwargs.get('depth')
    if (not points): points = 3
    if (not tolerance): tolerance = 1e-9
    if (not depth): depth = 10    
    dt = [None]*iterations
    for i in range(iterations):
        start = timer()
        result, error = xp.integral(func, a, b, points=points, tolerance=tolerance, depth=depth)
        end = timer()
        dt[i] = (end - start)*10**6 # in us
    mean = statistics.mean(dt)
    sigma = statistics.stdev(dt, xbar=mean)
    print('Average time (in us): ',str(mean),'+-',str(sigma))
    return [result, error]
'''

def plot_integral_analysis(number, true_val, begin, func, a, b, **kwargs):
    points = kwargs.get('points')
    tolerance = kwargs.get('tolerance')
    depth = kwargs.get('depth')
    if (not points): points = 3
    if (not tolerance): tolerance = 1e-9
    if (not depth): depth = 10
    
    est_errors = [None]*(points-begin+1)
    tr_errors = [None]*(points-begin+1)
    for i in range(begin,points+1):
        result, error = xp.integral(func, a, b, points=i, tolerance=tolerance, depth=depth)
        est_errors[i-begin] = error
        tr_errors[i-begin] = abs(true_val - result)
    
    pts = list(range(begin,points+1))

    plt.figure(number)

    #ploting
    plt.plot(pts, est_errors, label = "Estimated Error")
    plt.plot(pts, tr_errors, label = "True Error")

    plt.xlabel('Points')
    #plt.xticks(pts)
    plt.ylabel('Error')

    #plt.tick_params(axis='x',which='major',labelsize=3)

    plt.title('Performance Analysis')
    plt.legend()

    plt.show()

############# Function 1 #############

print("###Integrating the first function...\n")

phi = (1 + xp.sqrt(5)) / 2
acot = lambda x: xp.pi/2 - xp.atan(x)
true_res = 4 * xp.pi * acot(xp.sqrt(phi))

a = 0
b = xp.pi

func = lambda x: xp.tan(x) / xp.tan(x/2) * xp.log((2*xp.cos(x)**2 + 2*xp.cos(x) + 1) / (2*xp.cos(x)**2 - 2*xp.cos(x) + 1))

[result, error] = xp.integral(func, a, b, points=15)

print("Result:", result,"\nEst. Error:", error,"\nTrue Error:", result - true_res)

plot_integral_analysis(1, true_res, 3, func, a, b, points=100)

############# Function 2 #############

print("\n###Integrating the second function...\n")

true_res = (xp.pi / 8) * xp.log(xp.pi**2 / 8)

a = 0
b = 1

func = lambda x: xp.atan( (xp.atanh(x) - xp.atan(x)) / (xp.pi + xp.atanh(x) - xp.atan(x)) ) / x

[result, error] = xp.integral(func, a, b, points=10, tolerance=1e-15, depth=50)

print("Result:", result,"\nEst. Error:", error,"\nTrue Error:", result - true_res)

plot_integral_analysis(2, true_res, 3, func, a, b, points=100, tolerance=1e-15, depth=50)

############# Function 3 #############

print("\n###Integrating the third function...\n")

true_res = xp.pi / 2

a = 0
b = xp.inf

func = lambda x: (x**8 - 4*x**6 + 9*x**4 - 5*x**2 + 1) / (x**12 - 10*x**10 + 37*x**8 - 42*x**6 + 26*x**4 - 8*x**2 + 1)

[result, error] = xp.integral(func, a, b, points=10, tolerance=1e-15)

print("Result:", result,"\nEst. Error:", error,"\nTrue Error:", result - true_res)

plot_integral_analysis(3, true_res, 4, func, a, b, points=100, tolerance=1e-15)