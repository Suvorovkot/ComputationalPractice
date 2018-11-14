import numpy as np
import networkx as nx
import plotly as py
import plotly.graph_objs as go

eps = 3
n1 = 7
lenFr1 = np.array([[1062., 875, 650, 253, 421, 500, 250],
                   [875., 1027, 655, 235, 511, 499, 260],
                   [650., 655, 1039, 295, 570, 494, 256],
                   [253., 235, 295, 365, 254, 186, 115],
                   [421., 511, 570, 254, 814, 486, 222],
                   [500., 499, 494, 186, 486, 626, 221],
                   [250., 260, 256, 115, 222, 221, 452]])

lenFr2 = np.array([[548., 437, 286, 305, 313],
                   [437., 1020, 286, 286, 578],
                   [286., 286, 300, 286, 260],
                   [305., 286, 286, 474, 262],
                   [313., 578, 260, 262, 785]])
n2 = 5

def divLenghts(n, lenFr):
    divLen = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            divLen[i][j] = round(lenFr[i][j] / lenFr[i][i], eps)
    return divLen


def weights(n, divLen):
    wts = np.zeros((n, n))
    for j in range(n):
        for k in range(j + 1, n):
            wt = max(divLen[j][k], divLen[k][j])
            if(wt >= 0.5):
                wts[j][k] = wt
    return wts

def printWeights(n, weights):
    for j in range(n):
        for k in range(j, n):
            if(weights[j][k]!=0):
                print(j+1, "-", k+1,":", weights[j][k])
print("Debug:")
w1 = weights(n1, divLenghts(n1, lenFr1))
printWeights(n1, w1)

print("Fixup:")
w2 = weights(n2, divLenghts(n2, lenFr2))
#print(w2)
printWeights(n2, w2)







