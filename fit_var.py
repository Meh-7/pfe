import numpy as np
def fit_var(data, p):
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
    #residus
    residuals = X[p:] - X_hat
    return residuals, beta