import math
from prettytable import PrettyTable

print('\n', "Algebraic interpolating. Interpolating Lagrange and Newton polynomial", '\n')

def f(x):
    return math.exp(-x)-x**2/2

print('\n', "Test task: m = 15, x = 0.6, n = 7", '\n')

print('\n', "Please, enter number of nods", '\n')
m = int(input())
a = 0; b = 1

#Nods initialising
nods = []
for i in range(0, m+1):
    nods.append(a + (b - a) * i / m)

#Creating table
Init = PrettyTable()
Init.field_names = ["i", "Xi", "f(Xi)"]
for i in range(0, m+1):
    Init.add_row([i, nods[i], f(nods[i])])
print(Init, 'n')

print('\n', "Please, enter x value", '\n')
x = float(input())
print('\n', "Please, enter degree of interpolating polynomial", '\n')
n = int(input())

#Interpolation
def sortByResidual(n):
    return x - n
nods.sort(key = sortByResidual)
print(nods, '\n')
