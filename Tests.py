import numpy
import copy

from scipy import linalg

import GaussianMethod
import IterationMethodsLS

A = numpy.array([[8.29381, 0.995516, -0.560617],
                [0.995516, 6.298198, 0.595772],
                [-0.560617, 0.595772, 4.997407]])

b = numpy.array([[1], [0], [0]])

k = 7

xs = numpy.zeros((1, len(A[:, 0])))[0, :]
print("A|b:")
print(numpy.hstack((A, b)), "\n")

x = GaussianMethod.gaussFunc(numpy.hstack((A, b)))
print("1) Solution x*:")
print(x, "\n")

H, g = IterationMethodsLS.matrTransform(A, b)
print("2) H:")
print(H, "\n")
print("g:")
print(g, "\n")
print("||H|| = ", linalg.norm(H, numpy.inf), "\n")

est = IterationMethodsLS.aprEst(xs, H, g, k)
print("3) ||x^(7)-x*|| <= ", est, "\n")

xk = IterationMethodsLS.simpIter(xs, H, g, k)
print("4) x^(7):", xk, "\n")
print("||x^(7)-x*|| = ", linalg.norm(xk - x, numpy.inf), "\n")

xj = IterationMethodsLS.simpIter(xs, H, g, k-1)
est = IterationMethodsLS.apostEst(xk, xj, H)
print("||x^(7)-x*|| <= ", est, "\n")

L = IterationMethodsLS.Lusternik(H, xk, xj)
print("Lusternik's refinement:")
print(L)
print("||L(x^(7))-x*|| = ", linalg.norm(L - x, numpy.inf), "\n")

print("5) Seidel method:")
S = IterationMethodsLS.Seidel(H, g, xs, k)
print(S)
print("||xS-x*|| = ", linalg.norm(S - x, numpy.inf), "\n")

print("7) Relax method:")
print(IterationMethodsLS.sor(H, g, xs, k))


