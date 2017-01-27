# imports
from Protein import *
from Amino_acid import *
from algorithms import *
from child_protein import *
from time import sleep


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
# wave_search(sequence36, [(15,15)], 10, 400)

self_avoiding_walk(sequence14, 10)


