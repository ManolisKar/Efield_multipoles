import numpy as np
from lmfit import Model

def func(x,par0,par1,par2):
    return par0 + par1*x + par2*x**2

model = Model(func)
params = model.make_params(par0=0, par1=1, par2=0.25)

i, j = 7, 12
xdata = np.array([0] + x[:, i, j])
ydata = np.array([0] + y[:, i, j])

# now fit model to data, get results
result = model.fit(params, ydata, x=xdata)

print(result.fit_report())
