from scipy.misc import derivative
from scipy.misc import factorial as fac
import scipy.integrate as integrate
import numpy as np
from prettytable import PrettyTable


print('\n', "Finding solution for Cauchy problem for first-order ODE", '\n')

def f(y):
    return (-3) * y + y ** 2

h = 0.01
N = 10
x0 = 0
y0 = [1., -2., 2., 6., -30., -42.]

def solution(x):
    return 3 / (np.exp(3*x + np.log(2)) + 1)

def Taylor(x, x0, n):
    p = y0[0]
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
for j in range(5):
    print("|y", j, "- y(x",j ,")| = ", abs(fprimes[j]-solution(nods[j])))

#Adding Adams
print('\n', "Adams method", '\n')
Init = PrettyTable()
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(5, N+3):
    fprimes.append(fprimes[i-1] + 1/720 * h * (1901*f(fprimes[i-1]) - 2774*f(fprimes[i-2]) + 2616*f(fprimes[i-3]) - 1274*f(fprimes[i-4]) + 251*f(fprimes[i-5])))
    Init.add_row([i-2, nods[i], fprimes[i]])
print(Init, '\n')
print("|y10 - y(x10)| = ", abs(fprimes[10]-solution(nods[10])))

#Reinitialise nods
nods = []
for k in range(N+1):
    nods.append(x0 + k * h)

#Euler method
print('\n', "Euler's method", '\n')
Init = PrettyTable()
fprimes = [y0[0]]
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(1, N+1):
    fprimes.append(fprimes[i-1] + h * f(fprimes[i-1]))
    Init.add_row([i, nods[i], fprimes[i]])
print(Init, '\n')
print("|y10 - y(x10)| = ", abs(fprimes[9]-solution(nods[9])))

#Improvement
print('\n', "Improved Euler's method", '\n')
Init = PrettyTable()
fprimes = [y0[0]]
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(1, N+1):
    fprimes.append(fprimes[i-1] + (h/2) * (f(fprimes[i-1]) + f(fprimes[i-1] + h*f(fprimes[i-1]))))
    Init.add_row([i, nods[i], fprimes[i]])
print(Init, '\n')
print("|y10 - y(x10)| = ", abs(fprimes[9]-solution(nods[9])))

#Improvement2
print('\n', "Improved Euler's method 2", '\n')
Init = PrettyTable()
fprimes = [y0[0]]
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(1, N+1):
    fprimes.append(fprimes[i-1] + h * f(fprimes[i-1] + (h/2)*f(fprimes[i-1])))
    Init.add_row([i, nods[i], fprimes[i]])
print(Init, '\n')
print("|y10 - y(x10)| = ", abs(fprimes[9]-solution(nods[9])))


#R-K method
print('\n', "Runge-Kutt method", '\n')
Init = PrettyTable()
fprimes = [y0[0]]
Init.field_names = ["i", "Xi", "y(Xi)"]
for i in range(1, N+1):
    k1 = h*f(fprimes[i-1])
    k2 = h*f(fprimes[i-1] + (k1/2))
    k3 = h*f(fprimes[i-1] + (k2/2))
    k4 = h*f(fprimes[i-1] + k3)
    fprimes.append(fprimes[i-1] + (1/6) * (k1 + 2*k2 + 2*k3 + k4))
    Init.add_row([i, nods[i], fprimes[i]])
print(Init, '\n')
print("|y10 - y(x10)| = ", abs(fprimes[9]-solution(nods[9])))