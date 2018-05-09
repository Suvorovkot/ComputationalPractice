#import math
import np
import warnings
from prettytable import PrettyTable

warnings.simplefilter('ignore', np.RankWarning)

print('\n', "Inverse interpolation. Finding the derivatives of a table-defined function", '\n')

def f(x):
    return np.exp(-x) - x**2/2
def derf(x):
    return - np.exp(-x) - x
def secDerf(x):
    return np.exp(-x) - 1

def approxDer1(nods, i, h):
    return (nods[i+1][1] - nods[i-1][1]) / (2 * h)

def approxDer2(nods, h):
    return ((-3) * nods[0][1] + 4 * nods[1][1] - nods[2][1]) / (2 * h)

def approxDer3(nods, i, h):
    return (3 * nods[i][1] - 4 * nods[i-1][1] + nods[i-2][1]) / (2 * h)

def approxDer4(nods, i, h):
    return (nods[i+1][1] - 2 * nods[i][1] + nods[i-1][1]) / (h ** 2)

def polLagrange(deg, x, nds):
    P = 0
    for j in range(deg):
        p1 = 1; p2 = 1
        for i in range(deg):
            if i != j:
                p1 *= (x - nds[i][0])
                p2 *= (nds[j][0] - nds[i][0])
        P += nds[j][1] * p1 / p2
    return P

def dividedResidue(args):
    if len(args) == 1:
        return args[0][0]
    else:
        return (dividedResidue(args[1:][:]) - dividedResidue(args[0:-1][:])) / (args[-1][1] - args[0][1])

def polNewton(coef, x, nds):
    v = coef[0]
    for i in range(1, len(nds)):
        for j in range(0, i):
            coef[i] *= (x - nds[j][1])
        v += coef[i]
    return v

print('\n', "Please, enter number of nods", '\n')
m = int(input())

print('\n', "Please, enter [a,b]", '\n')
seg = input().split(",")
A = float(seg[0])
B = float(seg[1])

#Nods initialising
nods = []
for i in range(0, m + 1):
    x = (A + ((B - A) / m) * i)
    nods.append((x, f(x)))

#Creating initial table
Init = PrettyTable()
Init.field_names = ["i", "Xi", "f(Xi)"]
for i in range(0, m+1):
    Init.add_row([i, nods[i][0], nods[i][1]])
print(Init, '\n')

#First way (function is monotonous and continuous in [a,b])
print('\n', "Please, enter F value", '\n')
F = float(input())

print("--------------")
print("First way of inverse interpolation (f is monotonous and continuous in [", A, B, "])")
print('\n', "Please, enter degree of interpolating polynomial (<=", m,")", '\n')
deg = int(input())

while(deg > m):
    print("No, you're wrong, you're very wrong!")
    print("Please, enter degree of interpolating polynomial (<=", m,")", '\n')
    deg = int(input())


nods.sort(key=lambda nd: abs(F - nd[1]))
#print("Nods (sorted close to F):", '\n', nods)
coefNewton = []
for i in range(1, deg+1):
    coefNewton.append(dividedResidue(nods[:i][:]))
Q = polNewton(coefNewton, F, nods[:deg][:])
print("X = Qn(F) = ", Q)
print("|F - f(X)| = ", abs(F - f(Q)))
print("--------------", '\n')

print("--------------")
print("Second way of inverse interpolation")
#Second way (we don't know much about monotone)
nods.sort(key=lambda nd: nd[0])
print("Please, enter accuracy for finding roots", '\n')
e = float(input())

#RootSep
a = A; b = B
n = 50.0
h = (B - A) / n
x1 = a; x2 = x1 + h
k = 0
ints = []
while(x2 <= b):
    y1 = polLagrange(deg, x1, nods) - F
    y2 = polLagrange(deg, x2, nods) - F
    if (y1 * y2 <= 0):
        k += 1
        ints.append((x1, x2))
    x1 = x2
    x2 += h
#print('\n',"Number of intervals:", k)
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
    while abs(polLagrange(deg, x, nods) - F) >= e:
        if((polLagrange(deg, a, nods) - F) * (polLagrange(deg, x, nods) - F) < 0):
            a, b = a, x
            x = (a + b) / 2
        elif((polLagrange(deg, a, nods) - F) == 0):
            x = a
        else:
            a, b = x, b
            x = (a + b) / 2
        #count += 1
    #print("x",i,' = ',x, sep="")
    print("X = Qn(F) = ", x)
    print("|F - f(X)| = ", abs(F - f(x)))
    #print('\n',"With",count,"iterations.",'\n')
print("--------------", '\n')

print("--------------")
print("Finding derrivatives")
#print('\n', "Please, enter number of steps", '\n')
#s = int(input())

h = (B - A) / m

#Creating derrrivatives table
Init = PrettyTable()
Init.field_names = ["i", "Xi", "f(Xi)", "f'_a(Xi)", "|f'(Xi)-f'_a(Xi)|", "f''_a(Xi)", "|f''(Xi)-f''_a(Xi)|"]
Init.add_row([0, nods[0][0], nods[0][1], approxDer2(nods, h), abs(derf(nods[0][0])-approxDer2(nods, h)),
              "-----", "-----"])
for i in range(1, m):
    Init.add_row([i, nods[i][0], nods[i][1], approxDer1(nods, i, h), abs(derf(nods[i][0])-approxDer1(nods, i, h)),
                  approxDer4(nods, i, h), abs(secDerf(nods[i][0])-approxDer4(nods, i, h))])
Init.add_row([m, nods[m][0], nods[m][1], approxDer3(nods, m, h), abs(derf(nods[m][0])-approxDer3(nods, m, h)),
              "-----", "-----"])
print(Init, '\n')
