# imports
from Protein import *
from Amino_acid import *
from algorithms import *
from child_protein import *
from time import sleep
from math import *




# benchmark seq
sequence8 = ["H","H","P","H","H","H","P","H"]
sequence14 = ["H","H","P","H","H","H","P","H","P","H","H","H","P","H"]
sequence20 = ["H","P","H","P","P","H","H","P","H","P","P","H","P","H","H","P","P","H","P","H"] 
sequence36 = ["P","P","P","H","H","P","P","H","H","P","P","P","P","P","H","H","H","H","H","H","H","P","P","H","H","P","P","P","P","H","H","P","P","H","P","P"]




# locations = [(1,1,0),(0,1,0),(0,1,1),(1,1,1),(2,1,1),(2,1,2),(2,1,3),(2,1,4),(2,2,4),(2,2,3),(2,2,2)]

## SETTINGS TO REMEMBER

# complete doorrekening van seq14 = stability 6

#df(sequence8,[(15,15)])
#df(sequence14,[(15,15)])
# df(sequence8,[(15,15)])


# beam search solves seq 20 to max score of 9 (see paper Liu, Li, Yu)    BUT NOT ALWAYS?
#beam_search(sequence14, [(15,15)], 10, 1000)
#beam_search2(sequence36, [(15,15)], 2000)
# beam_search(sequence36, [(15,15)], 10, 1000)
#beam_search2(sequence20, [(15,15)], 2000)
# wave_search(sequence36, [(15,15)], 10, 400)


#self_avoiding_walk(sequence36, 10000)
# self_avoiding_walk(sequence14, 10)



############################ MOVE CHECKER
# protein = Protein(["H", "H", "H", "H", "H", "H", "H"],[(5,5),(5,6),(4,6),(4,5),(3,5),(2,5),(2,4)])
# protein.find_neighbours()
# protein.calculate_stability()
# protein.visualize()

# child = protein
# while True:
# 	options = [90, 180, 270]
# 	n = randint(2, 6)
# 	ang = randint(0,2)
# 	child.rotate(n, radians(options[ang]))
# 	child.find_neighbours()
# 	print child.validity_check()
	
# 	child.calculate_stability()

# 	if child.validity_check() == True:
# 		break

# child.visualize()

# options = [90, 180, 270]

# n = randint(2, 9)	
# ang = randint(0,2)
# protein.rotate(n,radians(options[ang]))

# protein.find_neighbours()
# protein.calculate_stability()

# print protein.validity_check()
# protein.visualize()
# protein.output()


######## HILLCLIMB

# # protein = Protein(sequence14,[(5,5),(5,6),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),(5,13),(5,14),(5,15),(5,16),(5,17),(5,18)])
# protein = Protein(sequence36,[(5,5),(5,6),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),(5,13),(5,14),(5,15),(5,16),(5,17),(5,18),(5,19),(5,20),(5,21),(5,22),(5,23),(5,24),(5,25),(5,26),(5,27),(5,28),(5,29),(5,30),(5,31),(5,32),(5,33),(5,34),(5,35),(5,36),(5,37),(5,38),(5,39),(5,40)])
# # protein = initialize(sequence14)
# protein.find_neighbours()
# protein.calculate_stability()

# hillclimbed = hillclimb(protein, 2000)
# protein = Protein(sequence14,[(5,5),(5,6),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),(5,13),(5,14),(5,15),(5,16),(5,17),(5,18)])
# protein = Protein(sequence20,[(5,5),(5,6),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),(5,13),(5,14),(5,15),(5,16),(5,17),(5,18),(5,19),(5,20),(5,21),(5,22),(5,23),(5,24)])
protein = Protein(sequence36,[(5,5),(5,6),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),(5,13),(5,14),(5,15),(5,16),(5,17),(5,18),(5,19),(5,20),(5,21),(5,22),(5,23),(5,24),(5,25),(5,26),(5,27),(5,28),(5,29),(5,30),(5,31),(5,32),(5,33),(5,34),(5,35),(5,36),(5,37),(5,38),(5,39),(5,40)])
# # protein = initialize(sequence14)

hillclimbed = simulated_annealing(protein, 2500)
# hillclimbed.find_neighbours()
# hillclimbed.calculate_stability()
# # hillclimbed.output()
hillclimbed.visualize()




