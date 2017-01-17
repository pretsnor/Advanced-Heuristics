import random
from Protein import *

def initialize():
	"""
	Generates a pseudorandom self avoiding walk of the input sequence
	"""

	sequence = ["H", "P", "H", "H", "P","H", "P", "H", "H", "P", "H"]
	locations = [(5,5,5)]

	for i in range(len(sequence)):
		while True:
			n = random.randint(0, 6)
			if n == 1:
				next = (locations[i][0],locations[i][1] - 1,locations[i][2])
			elif n == 2:
				next = (locations[i][0],locations[i][1] + 1,locations[i][2])
			elif n == 3:
				next = (locations[i][0],locations[i][1],locations[i][2] - 1)
			elif n == 4:
				next = (locations[i][0],locations[i][1],locations[i][2] + 1)
			elif n == 5:
				next = (locations[i][0] - 1,locations[i][1],locations[i][2])
			else:
				next = (locations[i][0] + 1,locations[i][1],locations[i][2])
			
			# avoids self intersections
			if next not in locations:
				break

		locations.append(next)

	# generate and visualize	
	protein = Protein(sequence,locations)
	protein.find_neighbours()
	protein.calculate_stability()

	protein.visualize()

