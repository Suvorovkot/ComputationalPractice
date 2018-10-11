import numpy
import copy

from scipy import linalg

import GaussianMethod




def matrTransform(A, b):
    l = len(A[:, 0])
    E = numpy.eye(l)
    D = numpy.zeros((l, l))
    for i in range(0, l):
        D[i][i] = A[i][i]
    H = E - numpy.dot(linalg.inv(D), A)
    G = numpy.dot(linalg.inv(D), b)
    g = numpy.zeros((1, l))[0, :]
    for i in range(0, l):
        g[i] = G[i][0]
    return H, g

def aprEst(xs, H, g, k):
    nH = linalg.norm(H, numpy.inf)
    return (nH**k) * linalg.norm(xs) + (nH**k / (1-nH)) * linalg.norm(g)

def apostEst(xk, xj, H):
    nH = linalg.norm(H, numpy.inf)
    return (nH / (1-nH)) * linalg.norm(xk - xj)

def simpIter(xs, H, g, k):
    def iterator(k):
        if k == 0:
            return xs
        else:
            return numpy.dot(H, iterator(k-1)) + g
    return iterator(k)

def Lusternik(H, xk, xj):
    rH = specRad(H)
    return xj + (1 / (1-rH)) * (xk - xj)

def specRad(H):
    eigVal = (linalg.eig(H)[0]).real
    rH = numpy.amax(eigVal)
    return rH

def Seidel(H, g, xs, k):
    l = len(H[:, 0])
    Hl = numpy.zeros((l, l))
    Hr = numpy.zeros((l, l))
    E = numpy.eye(l)
    for i in range(0, l):
        Hr[i][i:] = H[i][i:]
        Hl[i][:i] = H[i][:i]
    He = linalg.inv(E - Hl)
    def iterator(k):
        if k == 0:
            return xs
        else:
            return numpy.dot(numpy.dot(He, Hr), iterator(k-1)) + numpy.dot(He, g)
    return iterator(k)

def sor(H, g, xs, k):
    l = len(H[:, 0])
    Hl = numpy.zeros((l, l))
    Hr = numpy.zeros((l, l))
    for i in range(0, l):
        Hr[i][i+1:] = H[i][i+1:]
        Hl[i][:i] = H[i][:i]
    D = numpy.zeros((l, l))
    for i in range(0, l):
        D[i][i] = H[i][i]
    sr = specRad(H)
    #q = 1
    q = 2. / (1 + (1 - sr**2) * 0.5)
    def iterator(k):
        xk = numpy.zeros((1, l))[0, :]
        if k == 0:
            return xs
        else:
            xj = iterator(k-1)
            for i in range(0, l):
                s1 = 0
                for j in range(0, i):
                    s1 += (H[i][j] * xk[j])
                s2 = 0
                for j in range(i, l):
                    s2 += (H[i][j] * xj[j])
                xk[i] = xj[i] + q * (g[i] - xj[i] + s1 + s2)
            return xk
    return iterator(k)