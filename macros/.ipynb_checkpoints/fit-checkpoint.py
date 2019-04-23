import numpy as np
from ROOT import TCanvas, TGraph

#Some data
x = np.arange(10)
y = x**2

#Canvas to plot on and graph
can1 = TCanvas("can1","can1",800,800)
gr = TGraph(x.size, x.astype(np.double),y.astype(np.double))

#Draw canvas and graph
canvas.Draw()
graph.Draw()

