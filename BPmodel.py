import numpy as np

def BPstarmodel(A, B, omeg, n, delt):
    Y = np.zeros((2,n))
    lam = np.zeros(2)
    for i in range(1,n):
        lam = omeg + np.matmul(A,lam) + np.matmul(B,Y[:,i-1])
        yint = BPstar(lam[0], lam[1], delt)
        Y[0,i] = yint[0]
        Y[1,i] = yint[1]
    return Y