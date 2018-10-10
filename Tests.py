import numpy
import copy

from scipy import linalg
import unittest

import GaussianMethod
import IterationMethodsLS
import JacobiEigenvalueMethod

class TestCompMethods(unittest.TestCase):
    t1A = numpy.array([[8.29381, 0.995516, -0.560617],
                    [0.995516, 6.298198, 0.595772],
                    [-0.560617, 0.595772, 4.997407]])

    t2A = numpy.array([[4.0, 3.0, 2.0],
                     [3.0, 5.0, 1],
                     [2.0, 1.0, 7.0]])



    tb = numpy.array([[1], [0], [0]])


    def testIterMethods():
        A = TestCompMethods.t1A
        b = TestCompMethods.tb

        print("A|b:")
        print(numpy.hstack((A, b)), "\n")
        k = 7
        xs = numpy.zeros((1, len(A[:, 0])))[0, :]


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
        R = IterationMethodsLS.sor(H, g, xs, k)
        print(R)
        print("||xS-x*|| = ", linalg.norm(R - x, numpy.inf), "\n")


    def testJacobi():
        A = TestCompMethods.t1A
        print("A:")
        print(A, '\n')
        eigvalJ, eigvecJ = JacobiEigenvalueMethod.Jacobi(A, 1.0e-09)
        eigval, eigvec = linalg.eig(A)
        print("Eigenvalues:")

        print(eigvalJ, "\n")
        print("||e-eJ*|| = ", linalg.norm(numpy.sort(eigvalJ) - numpy.sort(eigval), numpy.inf), "\n")
        #print("||e-eJ*|| = ", eigvalJ - numpy.sort(eigval), "\n")
        print("Eigenvectors:")
        print(eigvecJ, "\n")
        L = numpy.identity(len(A))
        for i in range(len(L)):
            L[i][i] = eigvalJ[i]
        print("||R = AX - LX|| = ", linalg.norm(numpy.dot(A, eigvec) - numpy.dot(L, eigvec), numpy.inf), "\n")
        #print("|R = AX - LX| = ", numpy.abs(numpy.dot(A, eigvec) - numpy.dot(L, eigvec)), "\n")


if __name__ == '__main__':
    # TestCompMethods.testIterMethods()
     TestCompMethods.testJacobi()