import os
import matplotlib.pyplot as plt

DIR_INPUT = "input"
DIR_OUTPUT = "output"

for dirs in [DIR_INPUT, DIR_OUTPUT]:
    if not os.path.exists(dirs):
        os.makedirs(dirs)

#color = plt.cm.get_cmap('Classic', 10) #unnecessary unless you want different colors

COLOR_SCHEME = {'c' : 'k',
           'koff' : 'r',
           'simple_fisher' : 'b',
           'numerical_fisher_sp' : 'r',
           'numerical_fisher' : 'c',
           'heuristic' : 'y'}
