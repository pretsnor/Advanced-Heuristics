class Bond(object):
	"""
		A Bond object represents a chemical bond between two amino acids. 
	"""

	def __init__(self, aa1, aa2):
		self.bond = (aa1.location,aa2.location)

		self.value = 0
		self.type = aa1.acid + aa2.acid
		


	def determine_value(self):
		""" 	
		Determines the value of a bond 
		"""

		amino_matrix = {'HH': 1, 'HP': 0, 'PP': 0, 'PH': 0}
		self.value = amino_matrix[self.type]


	def __eq__(self, other):
		return self.bond == other.bond
		