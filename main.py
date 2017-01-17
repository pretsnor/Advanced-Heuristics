# imports
from Protein import *
from Amino_acid import *
	

sequence = ["H", "P", "H", "H", "P","H", "P", "H", "H", "P", "H"]
locations = [(1,1,0),(0,1,0),(0,1,1),(1,1,1),(2,1,1),(2,1,2),(2,1,3),(2,1,4),(2,2,4),(2,2,3),(2,2,2)]

test = Protein(sequence, locations)

test.find_neighbours()
test.calculate_stability()

test.visualize()

