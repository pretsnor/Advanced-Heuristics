class Bond(object):
	"""
		A Bond object represents a chemical bond between two amino acids. 
	"""

	def __init__(self, aa1, aa2):
		self.bond = (aa1.location,aa2.location)

	def __eq__(self, other):
		return self.bond == other.bond
		