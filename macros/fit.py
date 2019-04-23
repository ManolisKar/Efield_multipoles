import numpy as np
from ROOT import TCanvas, TGraph

def f(x,par):
    return par[0] + par[1]*x + par[2]*x**2

from scipy import optimize

minimum = optimize.fmin(f, 1)
print minimum