import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

def obvEuler(A, y0, h, k):
    W = np.identity(2) + h*A
    yk = y0
    Yi = np.array([[y0], [], [], []])
    for k in range(k):
        Yi[k] = yk
        yk = np.dot(W, yk)
    return Yi

def Adams(A, y0, h):
    W1 = linalg.inv(np.identity(2)-(5 * h / 12)*A)
    W2 = (np.identity(2) + (2 * h / 3) * A)
    Y1 = np.dot(np.identity(2)+h*A, y0)
    Y2 = np.dot(np.dot(W1, W2), Y1) - np.dot(np.dot(W1, (h / 12) * A), y0)
    return [Y1, Y2]

