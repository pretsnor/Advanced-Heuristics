class Amino_acid(object):
	"""
		Class amino acid represent amino acids that form proteins. Types possible: H and P (so far)
	"""

	def __init__(self, index, acid, location):
		self.acid = acid
		self.index = index
		self.location = location	

	# def find_neighbours(self):