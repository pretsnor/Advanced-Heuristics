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
#df(sequence14,[(15,15)])

# beam search solves seq 20 to max score of 9 (see paper Liu, Li, Yu)    BUT NOT ALWAYS?
#beam_search(sequence14, [(15,15)], 10, 1000)
#beam_search2(sequence36, [(15,15)], 2000)
# beam_search(sequence36, [(15,15)], 10, 1000)
# beam_search2(sequence20, [(15,15)], 2000)
# wave_search(sequence36, [(15,15)], 10, 400)

self_avoiding_walk(sequence14, 10)



############################ MOVE CHECKER
protein = Protein(["H", "H", "P", "H", "H", "H", "P", "P", "H", "P"], [(5,5),(5,6),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),(5,13),(5,14)])
protein.visualize()

options = [90, 180, 270]

for i in range(0, 5):
	n = random.randint(2, 9)	
	ang = random.randint(0,2)
	protein.rotate(n,radians(options[ang]))
	# check validity enzo!
protein.visualize()
protein.output()




