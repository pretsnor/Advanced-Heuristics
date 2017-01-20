# imports
import copy
from Protein import *

class Protein_child(Protein):
	"""
	Upon initialization, receives a deepcopy of parent characteristics.
	"""
	def __init__(self, parent):
		# get parent characteristics
		self.seq = copy.deepcopy(parent.seq)
		self.length = len(self.seq)
		self.locations = copy.deepcopy(parent.locations)
		self.x = copy.deepcopy(parent.x)
		self.y = copy.deepcopy(parent.y)

		self.stability = 0

		self.total_bonds = []
		self.covalent_bonds = []
		self.other_bonds = []

		# generate list (of tuples (of tuples)) containing all covalent bonds
		for i in range(self.length - 1):
			self.covalent_bonds.append(Bond(self.seq[i], self.seq[i + 1]))

	def add_amino(self, amino):
		self.seq.append(amino)
		self.length += 1
		self.locations.append(amino.location)
		self.x.append(amino.location[0])
		self.y.append(amino.location[1])