import numpy as np
import matplotlib.pyplot as plt
import pdb

x = np.linspace(-100, 100, 250)
y = np.linspace(-100, 100, 250)
xv, yv = np.meshgrid(x, y)
xv = xv.flatten()
yv = yv.flatten()
xv = np.expand_dims(xv, axis=1)
yv = np.expand_dims(yv, axis=1)

p = np.concatenate((xv, yv), axis = 1)

v0 = 10
p0 = np.asarray([30, 35])
v1 = 15
p1 = np.asarray([40, 60])

color = np.zeros((xv.shape[0], 3))
t0 = np.sqrt(np.sum((p - p0)**2, axis=1))/v0 + 1
t1 = np.sqrt(np.sum((p - p1)**2, axis=1))/v1

color[t0 >= t1, :] = [0, 1, 0];
color[t0 < t1, :] = [0, 0, 1];
plt.figure()
plt.scatter(xv, yv, c=color, marker='+')
plt.scatter(p0[0], p0[1], c=[1, 0, 0])
plt.scatter(p1[0], p1[1], c=[1, 0, 0])
plt.axis('equal')
plt.show()