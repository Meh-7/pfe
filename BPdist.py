import numpy as np

def BPstar(lam1, lam2, delt):
    c = 1 - np.e**(-1)
    val_lim = 50
    #initialisation de la premiere var poisson et var unif
    Y1 = np.random.poisson(lam=lam1, size=1)
    u = np.random.uniform(0,1)
    #fonction de masse
    P = np.zeros(val_lim+1)
    for i in range (0,val_lim):
        P[i] = (lam2**i/(np.math.factorial(i))) * np.e**(-lam2)*(1+ delt*float((np.e**(-Y1)-np.e**(-c*lam1))*(np.e**(-i)-np.e**(-c*lam2))))
    #fonction de mass cumulée
    Pcum = np.zeros(val_lim+1)
    Pcum[0] = P[0]
    for i in range(1,val_lim+1):
        Pcum[i] = P[i] + Pcum[i-1]
    
    found = False

    i = 0
    while (found == False) and (i<val_lim):
        if (u < Pcum[i]): 
            found = True
        else:
            i = i+1

    Y2 = np.array([i])
    Y = np.zeros((2,1))
    Y[0,0] = Y1
    Y[1,0] = Y2
    return Y  
    