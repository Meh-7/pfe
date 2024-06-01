import numpy as np
from fit_var import fit_var
from est_varma import estimate_varma
from myfunctions import center_data

def estimate_model(data, p): #p l'ordre du fit VAR
    cdata = center_data(data) #on centre les données

    resid, beta = fit_var(cdata, p) #on fit le var sur les données centrées et on tire les residus, beta optionnel
    omega, Theta, Phi = estimate_varma(data, p, resid) #on tire les regresseurs par les données originales, p et les residus

    #transformations pour A et B:
    A = -Theta
    B = Phi+Theta

    return omega, A, B
