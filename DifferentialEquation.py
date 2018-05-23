from scipy.misc import derivative
from scipy.misc import factorial as fac
import scipy.integrate as integrate
import numpy as np
from prettytable import PrettyTable


print('\n', "Finding solution for Cauchy problem for first-order ODE", '\n')

def f(y):
    return (-3) * y + y ** 2

h = 0.1
N = 10
x0 = 0
y0 = [1., -2., 2., 6., -30., -42.]

def solution(x):
    return 3 / (np.exp(3*x + np.log(2)) + 1)

def Taylor(x, x0, n):
    p = 1
    for i in range(1, n+1):
        p += ((y0[i]) / (fac(i))) * (x-x0)**i
    return p



#Creating table for accurate solution
print('\n', "Accurate solution", '\n')
nods = []
Init = PrettyTable()
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(-2, N+1):
    nods.append(x0 + i * h)
    Init.add_row([i, nods[i+2], solution(nods[i+2])])
print(Init, '\n')

#Creating table for Taylor
print('\n', "Taylor series method", '\n')
fprimes = []
Init = PrettyTable()
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(0, 5):
    fprimes.append(Taylor(nods[i], x0, 5))
    Init.add_row([(i-2), nods[i], fprimes[i]])
print(Init, '\n')

#Adding Adams
print('\n', "Adams method", '\n')
for i in range(5, N+1):
    fprimes.append(fprimes[i-1] + 1/720 * h * (1901*f(fprimes[i-1]) - 2774*f(fprimes[i-2]) + 2616*f(fprimes[i-3]) - 1274*f(fprimes[i-4]) + 251*f(fprimes[i-5])))
    Init.add_row([i-2, nods[i], fprimes[i]])
print(Init, '\n')

#Euler method
print('\n', "Euler's method", '\n')
Init = PrettyTable()
fprimes = [y0[0]]
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(0, N):
    fprimes.append(fprimes[i] + h * f(fprimes[i]))
    Init.add_row([i, nods[i], fprimes[i]])
print(Init, '\n')

#Improvement
print('\n', "Improved Euler's method", '\n')
Init = PrettyTable()
fprimes = [y0[0]]
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(0, N):
    fprimes.append(fprimes[i] + (h/2) * (f(fprimes[i]) + f(fprimes[i] + h*f(fprimes[i]))))
    Init.add_row([i, nods[i], fprimes[i]])
print(Init, '\n')

#R-K method
print('\n', "Runge-Kutt method", '\n')
Init = PrettyTable()
fprimes = [y0[0]]
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(0, N):
    k1 = h*f(fprimes[i])
    k2 = h*f(fprimes[i]+(k1/2))
    k3 = h*f(fprimes[i]+(k2/2))
    k4 = h*f(fprimes[i]+(k3))
    fprimes.append(fprimes[i] + (1/6) * (k1 + 2*k2 + 2*k3 + k4))
    Init.add_row([i, nods[i], fprimes[i]])
print(Init, '\n')
