import random
from Protein import *

def initialize():
	"""
	Generates a pseudorandom self avoiding walk of the input sequence
	"""

	sequence = ["H", "P", "H", "H", "P","H", "P", "H", "H", "P", "H","H"]
	locations = [(5,5)]

	for i in range(len(sequence)):
		while True:
			n = random.randint(0, 4)
			if n == 1:
				next = (locations[i][0],locations[i][1] - 1)
			elif n == 2:
				next = (locations[i][0],locations[i][1] + 1)
			elif n == 3:
				next = (locations[i][0] - 1,locations[i][1])
			else:
				next = (locations[i][0] + 1,locations[i][1])
			
			# avoids self intersections
			if next not in locations:
				break

		locations.append(next)

	# generate and visualize	
	protein = Protein(sequence,locations)
	
	return protein


	

