import numpy as np
from fit_var import fit_var
from est_varma import estimate_varma
def estimate_model(data, p): #p l'ordre du fit VAR
    resid, beta = fit_var(data, p) #on fit le var sur les données centrées et on tire les residus, beta optionnel
    omega, Theta, Phi = estimate_varma(data, p, resid) #on tire les regresseurs par les données originales, p et les residus

    #transformations pour A et B:
    A = -Theta
    B = Phi+Theta

    return omega, A, B
