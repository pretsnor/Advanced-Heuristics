# import
from Amino_acid import *
from Bond import *
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
		self.stability = 0

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

		self.total_bonds = []
		self.covalent_bonds = []
		self.other_bonds = []

		# generate list (of tuples (of tuples)) containing all covalent bonds
		for i in range(self.length - 1):
			self.covalent_bonds.append(Bond(self.seq[i], self.seq[i + 1]))


	def find_neighbours(self): 
		"""
		Finds neighbouring amino acid couples in the protein.
		"""	

		# find all possible bonds
		for i in range(self.length):
			for j in range(i, self.length):
				dx = abs(self.seq[i].location[0] - self.seq[j].location[0])
				dy = abs(self.seq[i].location[1] - self.seq[j].location[1])
				dz = abs(self.seq[i].location[2] - self.seq[j].location[2])
		
				if ((dx + dy + dz) == 1):
					self.total_bonds.append(Bond(self.seq[i], self.seq[j]))

		# subtract covalent bonds
		self.other_bonds = [x for x in self.total_bonds if x not in self.covalent_bonds]

	def calculate_stability(self):
		"""
		Calculates the stability of the protein in the current folding, based on neighbour interactions of amino acids
		"""
		for i in range(len(self.other_bonds)):
			self.other_bonds[i].determine_value()
			self.stability = self.stability + self.other_bonds[i].value

		print "stability of this fold: ", self.stability

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
		ax.set_xlabel('X')
		ax.set_ylabel('Y')
		ax.set_zlabel('Z')

		# plot amino acids
		ax.scatter(x, y, z, marker='o', color = colors)
		# plot covalent bonds
		ax.plot_wireframe(x, y, z, linestyle="-", color="#ff0000")
		
		# plot other bonds
		for i in range(len(self.other_bonds)):
			x = [self.other_bonds[i].bond[0][0], self.other_bonds[i].bond[1][0]]
			y = [self.other_bonds[i].bond[0][1], self.other_bonds[i].bond[1][1]]
			z = [self.other_bonds[i].bond[0][2], self.other_bonds[i].bond[1][2]]

			ax.plot_wireframe(x, y, z, linestyle=":", color="#000000")

		fig.suptitle('protein stability = %s'%(self.stability) , fontsize=14, fontweight='bold')


		plt.show()
 	
		
	