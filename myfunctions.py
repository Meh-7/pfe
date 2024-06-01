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

def BPstarmodel(A, B, omeg, n, delt):
    Y = np.zeros((2,n))
    lam = np.zeros(2)
    for i in range(1,n):
        lam = omeg + np.matmul(A,lam) + np.matmul(B,Y[:,i-1])
        yint = BPstar(lam[0], lam[1], delt)
        Y[0,i] = yint[0]
        Y[1,i] = yint[1]
    return Y


def fit_VAR(data, p):
    n, k = data.shape
    X = center_data(data)

    #matrice des regresseurs
    Xprime = np.zeros((n - p, k * p))
    X_target = np.zeros((n - p, k))

    for t in range(p, n):
        Xprime_row = []
        for lag in range(1, p + 1):
            Xprime_row.extend(X[t - lag])
        Xprime[t - p, :] = Xprime_row
        X_target[t - p] = X[t]
    
    #coefficients OLS
    beta = np.linalg.inv(Xprime.T @ Xprime) @ Xprime.T @ X_target
    #Prediction
    X_hat = Xprime @ beta
    residuals = X[p:] - X_hat
    return residuals, beta


def estimate_VAR_ARMA(data, p, residuals): #p = ordre du VAR sur les donnees centrees, necessaire pour les retards et valeurs manquantes des residus
    T, N = data.shape
    assert N == 2, "This implementation is for bivariate VAR only."
    print(f"T = {T}")
    n = T - p  # Number of rows after adjusting for lags
    Xprime = np.zeros((n-1, 5)) #n lignes pour les pbs, 3 = 1 cst + 2 termes retard X_t-1 + 2 terme retard et-1
    Xprime[:, 0] = 1  # colonne de la constante 
    print(f"n = {n}")
    for t in range(0, n-1):
        Xprime[t, 1:3] = data[p+t, :] #pour t= 0 le terme le plus lointain dans nos données (celui de retard 1) est le (p+1)eme element de nos données, et le premier qui a un terme de retard correspondant
        Xprime[t, 3:5] = residuals[t, :]
    print(Xprime)
    X1 = data[p+1:, 0]
    X2 = data[p+1:, 1]
    Xs = np.zeros((n-1,2))
    Xs[:,0] = X1
    Xs[:,1] = X2

    #Estimate coefficients using OLS
    beta = np.linalg.inv(Xprime.T @ Xprime) @ Xprime.T @ Xs

    Phi = np.zeros((2, 2))
    Theta = np.zeros((2, 2))
    omega = [beta[0, 0], beta[0, 1]]
    Theta[0, :] = beta[3,0], beta[4,0]
    Theta[1, :] = beta[3,1], beta[4,1]
    Phi[0, :] = beta[1,0], beta[1,1]
    Phi[1, :] = beta[2,0], beta[2,1]
    return omega, Theta, Phi