# import

from Amino_acid import *
from pprint import pprint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


class Protein(object):
	"""
		Class protein represents a protein. It consists of a sequence of amino acid objects.
	"""

	def __init__(self, sequence, locations):
		self.seq = []
		self.length = len(sequence)

		self.x = []
		self.y = []
		self.z = []

		# generate protein chain
		for i in range(self.length):
			# make amino acid
			self.seq.append(Amino_acid(i,sequence[i], locations[i]))

			# get locations
			self.x.append(locations[i][0])
			self.y.append(locations[i][1])
			self.z.append(locations[i][2])
	

	def __iter__(self):
		"""
			make protein object itarable
		"""
		return iter(self.seq)

	def visualize(self):
		"""
		Visualizes a protein in a 3d grid
		""" 

		# color coding
		color_code = {'H': 'red', 'P': 'blue'}
		colors = []

		# get colors of protein
		for i in range(self.length):
			colors.append(color_code[self.seq[i].acid])

		# coordinates
		x = self.x
		y = self.y
		z = self.z
		
		# plotting 	
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		# set limits
		ax.set_xlim(0,10)
		ax.set_ylim(0,10)
		ax.set_zlim(0,10)

		# axis labels
		ax.set_xlabel('X Label')
		ax.set_ylabel('Y Label')
		ax.set_zlabel('Z Label')

		# plot
		ax.scatter(x, y, z, marker='o', color = colors)
		ax.plot_wireframe(x, y, z, linestyle=":", color="#ff0000")

		plt.show()
 	
		
	# def calculate_stability(self):