import numpy
import copy

eps = 1e-16

def gaussFunc(a):

    a = numpy.array(a)

    len1 = len(a[:, 0])
    len2 = len(a[0, :])
    vectB = copy.deepcopy(a[:, len1])

    for i in range(len1):
        leadEl(a, i, len1)   #Chooses leading element in columns and swap lines

        alead = float(a[i][i])

        z = i
        while z < len2:
            a[i][z] = a[i][z] / alead
            z += 1

        j = i + 1
        while j < len1:
            b = a[j][i]
            z = i
            while z < len2:
                a[j][z] = a[j][z] - a[i][z] * b
                z += 1
            j += 1

    a = backTrace(a)

    #print("b - Ax:")
    #print(vectorN(c, a, len1, vectB))

    return a

def leadEl(a, i, len):
    max = abs(a[i][i])
    mi = i
    t1 = i
    while t1 < len:
        if abs(a[t1][i]) > max:
            max = abs(a[t1][i])
            mi = t1
        t1 += 1
    if abs(max) < eps:
        raise DetermExeption("Check matrix or use ")
    swapLines(a, i, mi)

def swapLines (a, i, j):
    if i != j:
        a[j], a[i] = copy.deepcopy(a[i]), copy.deepcopy(a[j])

class DetermExeption(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def backTrace(a):
    len1 = len(a[:, 0])
    len2 = len(a[0, :])

    a = numpy.array(a)
    i = len1 - 1
    while i > 0:
        j = i - 1
        while j >= 0:
            a[j][len1] = a[j][len1] - a[j][i] * a[i][len1]
            j -= 1
        i -= 1

    return a[:, len2 - 1]


def vectorN(c, a, rank, vectB):  # c = A, a - answer
    c = numpy.array(c)
    a = numpy.array(a)
    vectB = numpy.array(vectB)
    b = numpy.zeros((rank))
    i = 0
    while i < rank:
        j = 0
        while j < rank:
            b[i] += c[i][j]*a[j]
            j += 1
        i = i+1

    c = copy.deepcopy(b)

    for i in range(rank):
        c[i] = abs(c[i] - vectB[i])

    return c

def invMat(a):
    a = numpy.array(a)

    len1 = len(a[:, 0])

    E = numpy.eye(len1)
    a = numpy.hstack((a, E))
    len2 = len(a[0, :])
    for i in range(len1):
        #leadEl(a, i, len1)   #Chooses leading element in columns and swaps lines

        alead = float(a[i][i])

        z = i
        while z < len2:
            a[i][z] = a[i][z] / alead
            z += 1

        j = i + 1

        while j < len1:
            b = a[j][i]
            z = i
            while z < len2:
                a[j][z] = a[j][z] - a[i][z] * b
                z += 1
            j += 1

    for i in range(1, len1):
        k = i
        while k < len1:
            j = i - 1
            b = a[j][k]
            while j >= 0:
                z = k
                while z < len2:
                    a[j][z] -= a[k][z] * b
                    z += 1
                j -= 1
            k+=1
    return a


# a = numpy.array([[8.29381, 0.995516, -0.560617, 1],
#                 [0.995516, 6.298198, 0.595772, 0],
#                 [-0.560617, 0.595772, 4.997407, 0]])

# print("A|b:")
# print(a)
#
#
# d = invMat(a)
# print("\n")
# print("A^(-1):")
# for t in range(len(d[:, 0])):
#     print(d[t][len(d[0, :])-len(d[:, 0]):len(d[0, :])])
# print("\n")
#
# b = gaussFunc(a)
# print("\n")
# print("Answer:")
# print(b)
