import os
import matplotlib.pyplot as plt

DIR_INPUT = "input"
DIR_OUTPUT = "output"

for dirs in [DIR_INPUT, DIR_OUTPUT]:
    if not os.path.exists(dirs):
        os.makedirs(dirs)

color = plt.cm.get_cmap('Set1', 10)

COLOR_SCHEME = {'c' : 'k',
           'koff' : color(0),
           'simple_fisher' : color(1),
           'numerical_fisher_sp' : color(2),
           'numerical_fisher' : color(3),
           'heuristic' : color(4)}
