
import numpy
from numpy import linalg


def Jacobi(a, eps):
    def rotate(a, v, ik, jk):
        n = len(a)
        aDiv = a[ik][ik] - a[jk][jk]
        phi = 0.5 * numpy.arctan(2 * a[ik][jk] / aDiv)
        c = numpy.cos(phi)
        s = numpy.sin(phi)
        a[ik][jk] = 0.0
        a[jk][ik] = 0.0
        for i in range(n):
            if((i != ik)and(i != jk)):
                a[i][ik] = c * a[i][ik] + s * a[i][jk]
                a[ik][i] = a[i][ik]
                a[i][jk] = (-1) * s * a[i][ik] + c * a[i][jk]
                a[jk][i] = a[i][jk]
        a[ik][ik] = (c ** 2) * a[ik][ik] + 2 * c * s * a[ik][jk] + (s ** 2) * a[jk][jk]
        a[jk][jk] = (s ** 2) * a[ik][ik] - 2 * c * s * a[ik][jk] + (c ** 2) * a[jk][jk]
        for i in range(n):
            v[i][ik] = c * v[i][ik] + s * v[i][jk]
            v[i][jk] = (-1) * s * v[i][ik] + c * v[i][jk]


    n = len(a)
    v = numpy.identity(n) * 1.0
    for k in range(5):
        for i in range(n - 1):
            for j in range(i + 1, n):
                current = a[i][j]
                if numpy.abs(current) < eps: return numpy.diagonal(a), v
                rotate(a, v, i, j)
    print("Jacobi method did not converge")
    return numpy.diagonal(a), v
