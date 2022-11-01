import numpy as np
import matplotlib.pyplot as plt
import os

current_directory = os.path.dirname(os.path.realpath(__file__))

# read dataset
with open(os.path.join(current_directory, 'dataset.txt'), mode='r', encoding='utf-8') as f:
    lines = f.readlines()
    number_of_lines = len(lines)
    points_xyz = np.empty((number_of_lines, 3), dtype=np.float64)
    for line_no, line in enumerate(lines):
        for number_no, number in enumerate(line.split(' ')):
            points_xyz[line_no, number_no] = number.split('\n')[0]

centroid = points_xyz.mean(axis=0)
normal = np.array([0.330753, 0.543153, 0.771743], dtype=np.float64)

x_mean = round(centroid[0])
y_mean = round(centroid[1])

# create x,y
xx, yy = np.meshgrid(range(
    x_mean - 150, x_mean + 150), range(y_mean - 150, y_mean + 150))

# calculate corresponding z
d = -centroid.dot(normal)
z = (-normal[0] * xx - normal[1] * yy - d) * 1. / normal[2]

# Create the figure
fig = plt.figure()

# Add an axes
ax = fig.add_subplot(111, projection='3d')

# plot the surface
ax.plot_surface(xx, yy, z, alpha=0.3, color='aqua')


ax.scatter(points_xyz[:, 0], points_xyz[:, 1],
           points_xyz[:, 2],  color='green')

ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()
