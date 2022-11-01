import numpy as np
import matplotlib.pyplot as plt
import os

current_directory = os.path.dirname(os.path.realpath(__file__))

# define points at random:
p = np.random.randint(-100, 100, size=(50000, 3))

# A,B,C (c0) are arbitrary values; define D so plane intersects first point:
c0 = np.array([3, 5, 7])
D = -p[0].dot(c0)

# return all points in plane Ax + By + Cz + D = 0
plane_points = np.asarray(p[np.abs(p.dot(c0) + D) < 100], dtype=np.float64)
# print(plane_points)

# add bad outliers
noise = np.random.normal(-20, 20, size=(len(plane_points), 3))
for i in range(noise.shape[0]):
    for j in (2, 3, 5, 7):
        if (i % j) == 0:
            noise[i, :] = 0.0
            break
plane_points += noise

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter3D(plane_points[:, 0], plane_points[:, 1],
             plane_points[:, 2], 'gray', marker='.')
plt.show()

np.savetxt(os.path.join(current_directory, "dataset.txt"),
           plane_points, encoding='utf-8')
