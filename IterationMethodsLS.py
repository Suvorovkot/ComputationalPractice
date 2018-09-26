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
    def iterator (k):
        if k == 0:
            return xs
        else:
            return numpy.dot(H, iterator(k-1)) + g
    return iterator(k)