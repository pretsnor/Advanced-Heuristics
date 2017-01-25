# import
from Amino_acid import *
from Bond import *
from pprint import pprint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from math import *

class Protein(object):
	"""
		Class protein represents a protein. It consists of a sequence of amino acid objects.
	"""

	def __init__(self, sequence, locations):
		self.seq = []
		self.length = len(sequence)
		self.stability = 0
		self.locations = locations
		

		self.x = []
		self.y = []
		#self.z = []

		# generate protein chain
		for i in range(self.length):
			# make amino acid
			self.seq.append(Amino_acid(i,sequence[i], locations[i]))

			# get locations
			self.x.append(locations[i][0])
			self.y.append(locations[i][1])
			#self.z.append(locations[i][2])

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
				#dz = abs(self.seq[i].location[2] - self.seq[j].location[2])
		
				if ((dx + dy) == 1):
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

		#print "stability of this fold: ", self.stability

	def find_empty(self, location):
		"""
		Checks if there are empty spots around a location
		"""
		options = [(location[0] + 1, location[1]),
				   (location[0] - 1, location[1]),
				   (location[0], location[1] + 1),
				   (location[0], location[1] - 1)]

		# check if option is occupied by another amino acid
		empty_spots = [x for x in options if x not in self.locations]

		return empty_spots



	def rotate(self, n, angle):
	    """
	    Rotate all amino acids after amino acid n counterclockwise by a given angle around amino acid n.

	    The angle should be given in radians.

	    TODO: ZORGEN DAT IE DE KORTSTE TAK PAKT OM TE WISSELEN?
	    """

	    origin = self.seq[n].location
	    rest_length = (self.length - n)

	    ox, oy = origin
	    new_locs = []

	    # find all locations after
	    for i in range(n + 1, n + rest_length):
	    	px, py = self.seq[i].location
	    	qx = int(round(ox + cos(angle) * (px - ox) - sin(angle) * (py - oy)))
	    	qy = int(round(oy + sin(angle) * (px - ox) + cos(angle) * (py - oy)))
	    	# update all location sheit
	    	self.seq[i].location = (qx, qy)
	    	self.locations[i] = (qx, qy)
	    	self.x[i] = qx
	    	self.y[i] = qy

	def validity_check(self):
		if self.length != len(set(self.locations)):
			return False
		else:
			return True	


	def output(self):
		print "LOCATIONS: ", self.locations
		

	def __iter__(self):
		"""
		make protein object iterable
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
		#z = self.z
		
		# plotting 	
		fig = plt.figure()
		ax = fig.add_subplot(111)

		# set limits
		ax.set_xlim(0,30)
		ax.set_ylim(0,30)
		#ax.set_zlim(0,10)

		# axis labels
		ax.set_xlabel('X')
		ax.set_ylabel('Y')
		#ax.set_zlabel('Z')

		# plot amino acids
		ax.scatter(x, y, marker='o', color = colors)
		# plot grid
		#plt.grid(True,color='black')
		
		ax.plot(x, y, linestyle="-", color="#ff0000")
		
		# plot other bonds
		for i in range(len(self.other_bonds)):
			x = [self.other_bonds[i].bond[0][0], self.other_bonds[i].bond[1][0]]
			y = [self.other_bonds[i].bond[0][1], self.other_bonds[i].bond[1][1]]
			

			ax.plot(x, y, linestyle=":", color="#00ff00")

		fig.suptitle('protein stability = %s'%(self.stability) , fontsize=14, fontweight='bold')


		plt.show()
 	
		
	