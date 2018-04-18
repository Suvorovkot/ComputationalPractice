import math
from prettytable import PrettyTable

print('\n', "Algebraic interpolating. Interpolating Lagrange and Newton polynomial", '\n')

def f(x):
    return math.exp(-x)-x**2/2

def dividedResidue(f1, f2, x1, x2):
    return (f1 - f2) / (x1 - x2)

print('\n', "Test task: m = 15, x = 0.6, n = 7, a = 0, b = 1", '\n')

print('\n', "Please, enter number of nods", '\n')
#m = int(input())
m = 15
a = 0; b = 1

#Nods initialising
nods = []
for i in range(0, m+1):
    nods.append(a + (b - a) * i / m)

#Creating initial table
Init = PrettyTable()
Init.field_names = ["i", "Xi", "f(Xi)"]
for i in range(0, m+1):
    Init.add_row([i, nods[i], f(nods[i])])
print(Init, 'n')

print('\n', "Please, enter x value", '\n')
#x = float(input())
x = 0.6
print('\n', "Please, enter degree of interpolating polynomial (<=",m,")", '\n')
#deg = int(input())
deg = 7

#Newton interpolation
def sortByResidual(n):
    return abs(x - n)
nods.sort(key=sortByResidual)
Newton = []
#print(nods)

def dividedResidue(args):
    if len(args) == 1:
        return f(args[0])
    else:
        return (dividedResidue(args[1:]) - dividedResidue(args[0:-1])) / (args[-1] - args[0])

for i in range(1, deg+1):
    Newton.append(dividedResidue(nods[:i]))
print(Newton[:i])


