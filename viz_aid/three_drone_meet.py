import numpy as np
import matplotlib.pyplot as plt
import pdb

x = np.linspace(-200, 200, 250)
y = np.linspace(-200, 200, 250)
xv, yv = np.meshgrid(x, y)
xv = xv.flatten()
yv = yv.flatten()
xv = np.expand_dims(xv, axis=1)
yv = np.expand_dims(yv, axis=1)

p = np.concatenate((xv, yv), axis = 1)

q0 = np.asarray([20, 10])
qn = np.asarray([90, 60])

v0 = 10
p0 = np.asarray([10, 5])
v1 = 15
p1 = np.asarray([40, 60])
v1 = 25
p1 = np.asarray([70, 20])


color = np.zeros((xv.shape[0], 3))

plt.figure(1)
plt.scatter(xv, yv, c=color, marker='+')
plt.scatter(p0[0], p0[1], c=[1, 0, 0])
plt.scatter(p1[0], p1[1], c=[1, 0, 0])
plt.axis('equal')

plt.figure(2)
plt.scatter(xv, yv, c=color, marker='+')
plt.scatter(p0[0], p0[1], c=[1, 0, 0])
plt.scatter(p1[0], p1[1], c=[1, 0, 0])
plt.axis('equal')
plt.show()
plt.show()