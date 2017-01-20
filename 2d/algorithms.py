import random
from Protein import *
from child_protein import *

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

def df(sequence, start_position):
	"""
	Performs an exhaustive depth first search on a amino acid sequence, starting in start_position
	"""
	stack = []
	finished = []
	counter = 0	
	parent_counter = 0

	# start by generating 1 amino acid protein on the stack
	protein = Protein(sequence[0], start_position)
	stack.append(protein)

	for i in range(len(sequence) - 1):
		# pop parent from stack
		parent = stack.pop()
		parent_counter += 1

		# find empty spots next to last amino acid
		empty_spots = protein.find_empty(parent.seq[parent.length - 1].location)
		
		# generate childs with new amino acid on every possible spot
		for j in range(len(empty_spots)):
			protein = Protein_child(parent)
			amino = Amino_acid(i + 1, sequence[i + 1], empty_spots[j])
			protein.add_amino(amino)
			
			# count
			counter += 1

			# put new proteins on stack or in finished queue
			if protein.length < len(sequence):
				stack.append(protein)
			else:
				print "protein finished"
				finished.append(protein)
				






	

