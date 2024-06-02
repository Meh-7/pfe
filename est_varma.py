import numpy as np
def estimate_varma(data, p, residuals): #p = ordre du VAR sur les donnees centrees, necessaire pour les retards et valeurs manquantes des residus
    T, N = data.shape
    assert N == 2, "This implementation is for bivariate VAR only."
    #print(f"T = {T}")
    n = T - p  # Number of rows after adjusting for lags
    Xprime = np.zeros((n-1, 5)) #n lignes pour les pbs, 3 = 1 cst + 2 termes retard X_t-1 + 2 terme retard et-1
    Xprime[:, 0] = 1  # colonne de la constante 
    #print(f"n = {n}")
    for t in range(0, n-1):
        Xprime[t, 1:3] = data[p+t, :] #pour t= 0 le terme le plus lointain dans nos données (celui de retard 1) est le (p+1)eme element de nos données, et le premier qui a un terme de retard correspondant
        Xprime[t, 3:5] = residuals[t, :]
    
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