

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x =[1,2,2,2,3,6,7,8,8,8]
y =[6,6,7,8,8,8,8,8,8,8]
z =[0,0,0,0,0,0,0,1,1,2]
colors = ["red", "green", "blue", "red", "green", "blue", "red", "green", "blue", "red"]



ax.scatter(x, y, z, marker='o', color= colors)
ax.plot_wireframe(x, y, z, linestyle=":", color="#ff0000")

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')



plt.show()