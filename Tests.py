import numpy
import copy

import scipy
from scipy import linalg
import unittest

import GaussianMethod
import IterationMethodsLS
import JacobiEigenvalueMethod
import VariousEigenValuesMethods
import Adams
import DifferentialEquations2

class TestCompMethods(unittest.TestCase):
    t1A = numpy.array([[8.29381, 0.995516, -0.560617],
                    [0.995516, 6.298198, 0.595772],
                    [-0.560617, 0.595772, 4.997407]])

    t2A = numpy.array([[4.0, 3.0, 2.0],
                     [3.0, 5.0, 1],
                     [2.0, 1.0, 7.0]])

    t3A = numpy.array([[-125., 123.15],
                       [123.15, -123.]])

    t4A = numpy.array([[-125., 123.4],
                       [123.4, -123.]])

    t5A = numpy.array([[10.0, 2.0, 0., 0],
                     [3.0, 10.0, 4, 0],
                     [0.0, 1.0, 7.0, 5],
                     [0.0, 0.0, 3.0, 4]])

    tb1 = numpy.array([[1], [0], [0]])

    tb2 = numpy.array([3., 4, 5, 6])


    def testIterMethods(self):
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


    def testJacobi(self):
        A = TestCompMethods.t1A
        print("A:")
        print(A, '\n')
        eigvalJ, eigvecJ = JacobiEigenvalueMethod.Jacobi(A, 1.0e-09)
        eigval, eigvec = linalg.eig(A)
        print("Eigenvalues:")

        print(eigvalJ, "\n")
        # print("||e - eJ|| = ", linalg.norm(numpy.sort(eigvalJ) - numpy.sort(eigval), numpy.inf), "\n")
        print("|e - eJ| = ", numpy.sort(eigvalJ) - numpy.sort(eigval), "\n")
        print("Eigenvectors:")
        print(eigvecJ, "\n")
        L = numpy.identity(len(A))
        for i in range(len(L)):
            L[i][i] = eigvalJ[i]
        print("||R = AX - LX|| = ", linalg.norm(numpy.dot(A, eigvec) - numpy.dot(L, eigvec), numpy.inf), "\n")
        #print("|R = AX - LX| = ", numpy.abs(numpy.dot(A, eigvec) - numpy.dot(L, eigvec)), "\n")

    def testEigVal(self):
        eps = 0.001
        A = TestCompMethods.t1A
        eigval, eigvec = linalg.eig(A)
        print("A:")
        print(A, '\n')
        print("Python method:")
        print("lam1 = ", max(eigval), "\n")

        l1, x1, k = VariousEigenValuesMethods.powerMethod(A, eps)
        print("Power method:")
        print("lam1 = ", l1, " with eigenvector X = ", x1)
        print("||R = AX - lam1*X|| = ", linalg.norm(numpy.dot(A, x1) - l1*x1, numpy.inf), "with", k, "iterations \n")
        l1, x1, k = VariousEigenValuesMethods.scalProd(A, eps**2)
        print("Scalar composition method:")
        print("lam1 = ", l1, " with eigenvector X = ", x1)
        print("||R = AX - lam1*X|| = ", linalg.norm(numpy.dot(A, x1) - l1 * x1, numpy.inf), "with", k, "iterations \n")

        l, x = VariousEigenValuesMethods.specBound(A, eps, eigvec)
        print("Opposite bound of spectral radius is: ", l, "with eigenvector", x, "\n")

        l1, x1 = VariousEigenValuesMethods.Wielandt(A, eps)
        print("Wielandts method:")
        print("lam1 = ", l1, " with eigenvector X = ", x1)
        print("||R = AX - lam1*X|| = ", linalg.norm(numpy.dot(A, x1) - l1 * x1, numpy.inf), "\n")

    def testAdams(self):
        A = TestCompMethods.t4A
        # val, vec = linalg.eig(A)
        # print("!!!", val)
        h = 0.001
        k = 4
        y = numpy.array([1., 1.])
        print("Euler method:")
        for i in range(k):
            print("Y", i, " = ", Adams.obvEuler(A, y, h, k)[i])
        print("Interp Adams method:")
        for i in range(2):
            print("Y", i+1, " = ", Adams.Adams(A, y, h)[i])

    def testTDMA(self):
        A = TestCompMethods.t5A
        F = TestCompMethods.tb2
        print("Tridiagonal matrix method:")
        print(DifferentialEquations2.triDiagonalMatrixSolver(A, F))
        print("\n |Xm - X| = ", abs(DifferentialEquations2.triDiagonalMatrixSolver(A, F) - numpy.linalg.solve(A, F)))


if __name__ == '__main__':
    self = TestCompMethods()
    # TestCompMethods.testIterMethods(self)
    # TestCompMethods.testJacobi(self)
    # TestCompMethods.testEigVal(self)
    # TestCompMethods.testAdams(self)
    TestCompMethods.testTDMA(self)

