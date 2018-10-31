import numpy as np
from numpy import linalg


def powerMethod(A, eps):
    Yk = np.array([1., 0.0001, 0.0001])
    res = 1
    k = 0
    while res >= eps:
        k+=1
        Yk1 = np.dot(A, Yk)
        Yk1norm = np.linalg.norm(Yk1)
        Yk = Yk1 / Yk1norm
        l1 = Yk1[0]/Yk[0]
        res = linalg.norm(np.dot(A, Yk) - l1 * Yk, np.inf)
        if (k>=500):
            print("Probably does not converge")
            return 0, 0, k
    return l1, Yk, k

def scalProd(A, eps):
    Yk = np.array([1., 0.000001, 0.000001])
    res = 1
    k = 0
    while res >= eps:
        k += 1
        Yk1 = np.dot(A, Yk)
        Yk1norm = np.linalg.norm(Yk1)
        Yk = Yk1 / Yk1norm
        lam = np.dot(Yk1, Yk)/np.dot(Yk, Yk)
        res = linalg.norm(np.dot(A, Yk) - lam * Yk, np.inf)
        if (k>=500):
            print("Probably does not converge")
            return 0, 0, k
    return lam, Yk, k

def specBound(A, eps, evec):
    l, Y, k = scalProd(A, eps)
    B = A - l * np.identity(3)
    val, vec = linalg.eig(B)
    lB = min(val)
    evec[1] = evec[1] / linalg.norm(evec[1])
    return lB + l, evec[1]

def Wielandt(A, eps):
    lk = 8.
    k = 0
    res = 1
    while res >= eps:
        k += 1
        W = A - lk * np.identity(3)
        l, v, i = powerMethod(linalg.inv(W), eps)
        lk = 1/l + lk
        res = linalg.norm(np.dot(A, v) - lk * v, np.inf)
        if (k>=500):
            print("Probably does not converge")
            return 0, 0, k
    return lk, v
