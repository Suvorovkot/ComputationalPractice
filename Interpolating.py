import math
from prettytable import PrettyTable

print('\n', "Algebraic interpolation. Interpolating Lagrange and Newton polynomial", '\n')

def f(x):
    return math.exp(-x)-x**2/2

def dividedResidue(args):
    if len(args) == 1:
        return f(args[0])
    else:
        return (dividedResidue(args[1:]) - dividedResidue(args[0:-1])) / (args[-1] - args[0])

def polNewton(coef, x, nds):
    v = coef[0]
    for i in range(1, len(nds)):
        for j in range(0, i):
            coef[i] *= (x - nds[j])
        v += coef[i]
    return v

print('\n', "Test task: m = 15, x = 0.6, n = 7, a = 0, b = 1", '\n')

print('\n', "Please, enter number of nods", '\n')
m = int(input())
print('\n', "Please, enter [a,b]", '\n')
seg = input().split(",")
a = int(seg[0])
b = int(seg[1])

#Nods initialising
nods = []
for i in range(1, m+2):
    #nods.append(a + ((b - a) / 11) * (i - 1))
    nods.append(((b - a) / 2) * math.cos(((2 * i - 1) * math.pi) / (2 * (m + 1))) + ((b + a) / 2))

#Creating initial table
Init = PrettyTable()
Init.field_names = ["i", "Xi", "f(Xi)"]
for i in range(0, m+1):
    Init.add_row([i, nods[i], f(nods[i])])
print(Init, 'n')

print('\n', "Please, enter x value", '\n')
x = float(input())

print('\n', "Please, enter degree of interpolating polynomial (<=",m,")", '\n')
deg = int(input())
while(deg > m):
    print("Please, enter degree of interpolating polynomial (<=",m,")", '\n')


def sortByResidual(n):
    return abs(x - n)
nods.sort(key=sortByResidual)
print("Nods (sorted close to x):", '\n', nods[:deg])

#Newton interpolation
coefNewton = []
for i in range(1, deg+1):
    coefNewton.append(dividedResidue(nods[:i]))
P = polNewton(coefNewton, x, nods[:deg])
print("--------------")
print("Newton's: Pn(x) = ", P)
print("|f(x) - Pn(x)| = ", math.fabs(f(x)-P))
print("--------------",'\n')

#Lagrange interpolation
P = 0
fnods = []
for nd in nods[:deg]:
    fnods.append(f(nd))
for j in range(deg):
    p1 = 1; p2 = 1
    for i in range(deg):
        if i != j:
            p1 *= (x - nods[i])
            p2 *= (nods[j] - nods[i])
    P += fnods[j] * p1 / p2
print("--------------")
print("Lagrange's: Pn(x) = ", P)
print("|f(x) - Pn(x)| = ", math.fabs(f(x)-P))
print("--------------",'\n')
