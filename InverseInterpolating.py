#import math
import np
import warnings
from prettytable import PrettyTable

warnings.simplefilter('ignore', np.RankWarning)

print('\n', "Inverse interpolation. Finding the derivatives of a table-defined function", '\n')

def f(x):
    return np.exp(-x) - x**2/2
def derf(x):
    return -np.exp(-x) - x
def secDerf(x):
    return np.exp(-x) - 1

def approxDer1(y, h):
    return (f(y + h) - f(y - h)) / (2 * h)

def approxDer2(y, h):
    return (-3*f(y) + 4 * f(y + h) - f(y + 2 * h)) / (2 * h)

def approxDer3(y, h):
    return (3*f(y) - 4 * f(y + h) + f(y + 2 * h)) / (2 * h)

def approxDer4(y, h):
    return (f(y + h) - 2 * f(y) + f(y - h)) / (h ** 2)

def polLagrange(deg, x, nds, fnds):
    P = 0
    for j in range(deg):
        p1 = 1; p2 = 1
        for i in range(deg):
            if i != j:
                p1 *= (x - nds[i])
                p2 *= (nds[j] - nds[i])
        P += fnds[j] * p1 / p2
    return P

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

print('\n', "Please, enter number of nods", '\n')
m = int(input())

print('\n', "Please, enter [a,b]", '\n')
seg = input().split(",")
a = float(seg[0])
b = float(seg[1])

#Nods initialising
nods = []
for i in range(0, m + 1):
    nods.append(a + ((b - a) / m) * i)

#Creating initial table
Init = PrettyTable()
Init.field_names = ["i", "Xi", "f(Xi)"]
for i in range(0, m+1):
    Init.add_row([i, nods[i], f(nods[i])])
print(Init, '\n')

#First way (function is monotonous and continuous in [a,b])
print('\n', "Please, enter F value", '\n')
F = float(input())

print("--------------")
print("First way of inverse interpolation (f is monotonous and continuous in [", a, b, "])")
print('\n', "Please, enter degree of interpolating polynomial (<=",m,")", '\n')
deg = int(input())

while(deg > m):
    print("No, you're wrong, you're very wrong!")
    print("Please, enter degree of interpolating polynomial (<=",m,")", '\n')
    deg = int(input())

fnods = []
for nd in nods[:deg]:
    fnods.append(f(nd))
def sortByResidual(n):
    return abs(F - n)
fnods.sort(key=sortByResidual)
print("Nods (sorted close to F):", '\n', fnods[:deg])

coefNewton = []
for i in range(1, deg+1):
    coefNewton.append(dividedResidue(fnods[:i]))
Q = polNewton(coefNewton, F, fnods[:deg])
print("Newtons's: X = Qn(F) = ", Q)
print("|F - f(X)| = ", np.fabs(F-f(Q)))
print("--------------", '\n')

print("--------------")
print("Second way of inverse interpolation")
#Second way (we don't know much about monotone)
fnods = []
for nd in nods:
    fnods.append(f(nd))
print("Please, enter accuracy for finding roots", '\n')
e = float(input())

#RootSep
n = 50.0
h = (b-a)/n
x1 = a; x2 = x1 + h
k = 0
ints = []
while(x2 <= b):
    y1 = polLagrange(deg, x1, fnods, nods) - F
    y2 = polLagrange(deg, x2, fnods, nods) - F
    if (y1 * y2 <= 0):
        k += 1
        ints.append((x1,x2))
    x1 = x2
    x2 += h
print('\n',"Number of intervals:", k)
#for i in range(0, k):
#    print('[',ints[i][0],',',ints[i][1],']')
print('\n')

#Bisection
#print("----------------------",'\n',"Bisection method")
for i in range(0, k):
    a = ints[i][0]; b = ints[i][1]
    x = (a + b) / 2
    #print("For x", i, ' initial approximation is: ', x, sep="")
    #count = 1
    while np.fabs(polLagrange(deg, x, fnods, nods) - F) >= e:
        a, b = (a, x) if (polLagrange(deg, a, fnods, nods) - F) * (polLagrange(deg, x, fnods, nods) - F) < 0 else (x, b)
        x = (a + b) / 2
        #count += 1
    #print("x",i,' = ',x, sep="")
    print("Lagrange's: X = Qn(F) = ", x)
    print("|F - f(X)| = ", np.fabs(F - f(x)))
    #print('\n',"With",count,"iterations.",'\n')
print("--------------", '\n')

print("--------------")
print("Finding derrivatives")
h = (b - a) / m

#Creating derrrivatives table
Init = PrettyTable()
Init.field_names = ["i", "Xi", "f(Xi)", "f'_a(Xi)", "|f'(Xi)-f'_a(Xi)|", "f''_a(Xi)", "|f''(Xi)-f''_a(Xi)|"]
Init.add_row([0, nods[0], f(nods[0]), approxDer2(nods[0], h), abs(derf(nods[0])-approxDer2(nods[0], h)), approxDer4(nods[0], h), abs(secDerf(nods[0])-approxDer4(nods[0], h))])
for i in range(1, m):
    Init.add_row([i, nods[i], f(nods[i]), approxDer1(nods[i], h), abs(derf(nods[0])-approxDer2(nods[i], h)), approxDer4(nods[i], h), abs(secDerf(nods[i])-approxDer4(nods[i], h))])
Init.add_row([m, nods[m], f(nods[m]), approxDer3(nods[m], h), abs(derf(nods[m])-approxDer3(nods[m], h)), approxDer4(nods[m], h), abs(secDerf(nods[m])-approxDer4(nods[m], h))])
print(Init, '\n')
