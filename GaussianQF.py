import numpy
import scipy.integrate as integrate
import scipy.special.orthogonal
import warnings
import np
warnings.simplefilter('ignore', numpy.RankWarning)

print('\n', "Integration with Gaussian quadrature formulas", '\n')
print('\n', "Please, enter f(x)", '\n', "(Please use np prefix: np.exp, np.sin, np.cos and etc)", '\n')
fstr = input()
def check(x):
    return x**3
def f(x):
    return float(eval(fstr))

print('\n', "Please, enter w(x)", '\n', "(Please use np prefix: np.exp, np.sin, np.cos and etc)", '\n')
wstr = input()
def w(x):
    return float(eval(wstr))

def fw(x):
    return f(x)*w(x)

print('\n', "Please, enter [a,b]", '\n')
seg = input().split(",")
a = float(seg[0])
b = float(seg[1])

print('\n', "Please, enter m", '\n')
m = int(input())

print("--------------")
print("Moments of weight function:")
mts = []
for k in range(4):
    mts.append(integrate.quad(lambda x: w(x)*(x**k), a, b)[0])
    print("M" + str(k) + " =", mts[k])

print("--------------")
print("w2(x) = x^2 + a1*x + a2:")
ai = ((mts[0]*mts[3] - mts[1]*mts[2]) / (mts[1]**2 - mts[0]*mts[2]), (mts[2]**2 - mts[1]*mts[3]) / (mts[1]**2 - mts[0]*mts[2]))
print("a1 =", ai[0])
print("a2 =", ai[1])

print("--------------")
print("Roots for w2(x):")
cffs = [1, ai[0], ai[1]]
xi = np.roots(cffs)
print(xi)
print("x1 =", xi[0])
print("x2 =", xi[1])
print("--------------")
print("Coefficients for quadrature formula of Gaussian type:")
A1 = 1/(xi[0] - xi[1]) * (mts[1] - xi[1] * mts[0])
A2 = 1/(xi[1] - xi[0]) * (mts[1] - xi[0] * mts[0])
print("A1 =", A1)
print("A2 =", A2)
print("Checking: |M3 - A1(x1)**3 - A2(x2)**3| =", abs(mts[3]-A1*(xi[0])**3 - A2*(xi[1])**3))
print("--------------")
J = integrate.quad(lambda x: w(x)*f(x), a, b)
#print("J = ", J[0])
#print("--------------")
I = A1*f(xi[0]) + A2*f(xi[1])
print("I = ", I)
print("|J - I| = ", abs(I - J[0]))
print("--------------")
print("Gaussian quadrature formula:")
[nds, cfs] = scipy.special.roots_legendre(2)
print("Nods:", nds)
print("Coefficients:", cfs)
h = (b - a) / m
JG = 0
nds *= (h / 2)
for k in range (m):
    JG += (fw(nds[0] + a + k*h + (h / 2)) + fw(nds[1] + a + k*h + (h / 2)))
JG *= (h / 2)
print("I = ", JG)

print("|J - I| = ", abs(JG - J[0]))






