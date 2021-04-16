
import matplotlib.pyplot as plt
import numpy as np


@np.vectorize
def exponential_shelf(x, h=100, a=1e-5):
    return - h * (1. - np.exp(- 1. * a * x))


x = np.linspace(0, 1e6, 501)
y = np.linspace(-1e7, 1e7, 5)
xx, yy = np.meshgrid(x, y)
zz = exponential_shelf(xx)


plt.figure()
plt.plot(x, zz[0, :])
plt.show()
