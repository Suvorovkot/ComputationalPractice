import numpy


def triDiagonalMatrixSolver(A, F):
    n = len(F)
    a = numpy.zeros(n-1)
    b = numpy.zeros(n)
    c = numpy.zeros(n-1)
    d = F
    for i in range(n):
        for j in range(n):
            if(i==j):
                b[i] = A[i][j]
            elif(i+1 == j):
                c[i] = A[i][j]
            elif(i-1 == j):
                a[i-1] = A[i][j]
    ac, bc, cc, dc = map(numpy.array, (a, b, c, d))
    for it in range(1, n):
        mc = ac[it - 1] / bc[it - 1]
        bc[it] = bc[it] - mc * cc[it - 1]
        dc[it] = dc[it] - mc * dc[it - 1]

    xc = bc
    xc[-1] = dc[-1] / bc[-1]

    for il in range(n-2, -1, -1):
        xc[il] = (dc[il] - cc[il] * xc[il + 1]) / bc[il]

    return xc


def signing(x, y):
    if(y>=0):
        return abs(x)
    else:
        return abs(x)*(-1.)


def toTDM(A):    #Householder's method
    n = len(A)
    sx = numpy.zeros(n)
    for i in range(n-2):
        for k in range(i, n):
            sx[k] = A[n-1][i]*A[n-1][k]
        for j in range(n - 1, i + 1, -1):
            sx[i] += (A[j][i] * A[j][i])
        for k in range(i+1, n):
            for j in range(n-1, k, -1):
                sx[k] = sx[k] + (A[j][i]*A[j][k])
            for j in range(k-1, i+1, -1):
                sx[k] = sx[k] + (A[j][i]*A[k][j])

        alpha = numpy.emath.sqrt(sx[i])
        if (A[i+1][i]!=0):
            beta = 1./alpha
            for j in range(i+2, n):
                A[j][i] *= beta
            sx[i] = A[i+1][i]*beta + signing(1., A[i+1][i])
            A[i+1][i] = alpha
            G = 1. / abs(sx[i])
            sx2 = 0.
            for k in range(i+2,n):
                sx[k] = sx[k]*beta*G + signing(A[k][i+1],sx[i])
                sx2 += sx[k]*A[k][i]
            sx2 *= G
            for k in range(i+2, n):
                A[k][k] -= (2 * A[k][i] * sx[k] + (sx2 * A[k][i]) ** 2)
                for j in range(k+1,n):
                    A[j][k] -= (A[j][i] * sx[k] - A[k][i] * sx[I] + sx2 * A[j][i] * A[k][i])
        elif(alpha!=0):
            beta = 1. / alpha
            for j in range(i + 2, n):
                A[j][i] *= beta
            sx[i] = -1.
            A[i+1][i] = alpha
            G = 1.
            sx2 = 0.
            for k in range(i+2, n):
                sx[k] = sx[k] * beta * G + signing(A[k][i + 1], sx[i])
                sx2 += sx[k] * A[k][i]
            sx2 *= G
            for k in range(i+2, n):
                A[k][k] -= (2 * A[k][i] * sx[k] + (sx2 * A[k][i]) ** 2)
                for j in range(k+1,n):
                    A[j][k] -= (A[j][i] * sx[k] - A[k][i] * sx[I] + sx2 * A[j][i] * A[k][i])
        else:
            sx[i] = 1.
        Ans = numpy.zeros((n,n))
        for i in range(n):
            for j in range(n):
                if((i==j) or (i-1==j) or(i==j-1)):
                    Ans[i][j] = A[i][j]

    return Ans