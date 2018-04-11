import math

print('\n',"Algebraic interpolating. Interpolating Lagrange and Newton polynomial",'\n')

def f(x):
    return math.exp(-x)-x**2/2

m = 15; n = 7; a = 0; b = 1; x = 0.6

#Nods initialising
nods = []
for i in range(0,m+1):
    nods.append(a + (b - a) * i / m)
#print(nods)
