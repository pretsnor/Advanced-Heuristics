import random
from Protein import *
from child_protein import *
from heapq import *
from time import sleep
from Queue import *


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
	best = 0

	# start by generating 1 amino acid protein on the stack
	protein = Protein(sequence[0], start_position)
	stack.append(protein)

	while len(stack) > 0:
		# pop parent from stack
		parent = stack.pop()
		parent_counter += 1

		# find empty spots next to last amino acid
		empty_spots = protein.find_empty(parent.seq[parent.length - 1].location)
		
		# generate childs with new amino acid on every possible spot
		for j in range(len(empty_spots)):
			protein = Protein_child(parent)
			amino = Amino_acid(parent.length, sequence[parent.length], empty_spots[j])
			protein.add_amino(amino)
			
			# put new proteins on stack or in finished queue
			if protein.length < len(sequence):
				stack.append(protein)
			else:
				protein.find_neighbours()
				protein.calculate_stability()
				print "protein stab", protein.stability
				counter += 1
				print "protein number: ", counter
				# TO DO optimize: keep track of best one so far
				if protein.stability < best:
					heappush(finished,(protein.stability,protein))
					best = protein.stability
	
	# get best fold from prio queue
	best = heappop(finished)[1]
	best.visualize()	

def df_pruning(sequence, start_position):
	"""
	Performs an exhaustive depth first search on a amino acid sequence, starting in start_position
	"""
	stack = []
	finished = []
	counter = 0	
	parent_counter = 0
	best = 0

	# start by generating 1 amino acid protein on the stack
	protein = Protein(sequence[0], start_position)
	stack.append(protein)

	while len(stack) > 0:
		# pop parent from stack
		parent = stack.pop()
		parent_counter += 1

		# find empty spots next to last amino acid
		empty_spots = protein.find_empty(parent.seq[parent.length - 1].location)
		
		# generate childs with new amino acid on every possible spot
		for j in range(len(empty_spots)):
			protein = Protein_child(parent)
			amino = Amino_acid(parent.length, sequence[parent.length], empty_spots[j])
			protein.add_amino(amino)

			protein.find_neighbours()
			protein.calculate_stability()
			# put new proteins on stack or in finished queue
			if protein.length < 6:
				stack.append(protein)
			elif protein.length < len(sequence):
				if protein.stability > 1:
					print " NOT PRUNED"
					print "stab", protein.stability
					print "len", protein.length
					stack.append(protein)
				else:
					print "PRUNED "
					print "stab", protein.stability
					print "len", protein.length
			else:
				# TO DO optimize: keep track of best one so far
				if protein.stability > best:
					heappush(finished,((protein.stability * -1),protein))
					best = protein.stability
	
	# get best fold from prio queue
	best = heappop(finished)[1]
	best.visualize()	

def beam_search(sequence, start_position, k, n):
	"""
	Based on a breadth first algorithm. Every k steps, we take the best n proteins (resulting in a wave shape pruning).
	"""
	prioqueue = []


	finished = []
	counter = 0	
	parent_counter = 0
	best = 0

	# start by generating 1 amino acid protein on the stack
	protein = Protein(sequence[0], start_position)
	heappush(prioqueue,(0, protein))

	while len(prioqueue) > 0:

		# check for beam
		if (len(prioqueue) > k * n):
			print parent.length
			print "BEFORE: ", len(prioqueue)
			prioqueue = nsmallest(n,prioqueue)
			print "AFTER: ", len(prioqueue)

		# pop parent from prioqueue
		parent = heappop(prioqueue)[1]
		
		# find empty spots next to last amino acid
		empty_spots = parent.find_empty(parent.seq[parent.length - 1].location)
		
		# generate childs with new amino acid on every possible spot
		for j in range(len(empty_spots)):

			protein = Protein_child(parent)
			amino = Amino_acid(parent.length, sequence[parent.length], empty_spots[j])
			protein.add_amino(amino)

			protein.find_neighbours()
			protein.calculate_stability()
			# put new proteins on stack or in finished queue
		
			if protein.length < len(sequence):
				heappush(prioqueue,(protein.stability, protein))
			else:
				# TO DO optimize: keep track of best one so far
				if protein.stability <= best:
					counter += 1
					print "finished protein: ", counter
					heappush(finished,(protein.stability,protein))
					best = protein.stability
	
	# get best fold from prio queue
	best = heappop(finished)[1]

	best.output()
	best.visualize()	


def beam_search2(sequence, start_position, w):
	"""
	Based on a breadth first algorithm. Take w children every generation.

	this one is correcto

	todo: prune somehow
	"""

	total_queue = Queue()
	total_counter = 0

	temp_queue = []
	temp_counter = 0
	best = 0

	finished = []
	finished_counter = 0

	# make first anchestor protein
	ancestor = Protein(sequence[0], start_position)
	total_queue.put(ancestor)

	while (total_queue.qsize() > 0):		
		parent = total_queue.get()

		# check if finished
		if parent.length == len(sequence):
			if parent.stability <= best:
				heappush(finished,(parent.stability, finished_counter, parent))
				finished_counter += 1
				best = parent.stability
		else:

			# find empty spots next to last amino acid
			empty_spots = parent.find_empty(parent.seq[parent.length - 1].location)

			# make child with next amino acid on all empty spots
			for i in range(len(empty_spots)):
				child = Protein_child(parent)
				amino = Amino_acid(parent.length, sequence[parent.length], empty_spots[i])
				child.add_amino(amino)
				child.find_neighbours()
				child.calculate_stability()

				if (child.length == 10 and child.stability > -3):
					break
				else:
					print child.stability
					print child.length
					total_queue.put(child)
					temp_counter += 1

	best = heappop(finished)[2]
	best.visualize()
		








	

