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
