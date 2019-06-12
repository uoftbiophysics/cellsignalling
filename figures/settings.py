import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import seaborn as sns

DIR_INPUT = "input"
DIR_OUTPUT = "output"

for dirs in [DIR_INPUT, DIR_OUTPUT]:
    if not os.path.exists(dirs):
        os.makedirs(dirs)

# color = plt.cm.get_cmap('Classic', 10)  # unnecessary unless you want different colors
colour_palette = sns.color_palette("muted", 10)
#sns.palplot(colour_palette)
#plt.show()

COLOR_SCHEME = {'c' : 'k',
                'koff' : colour_palette[3],
                'simple_fisher' : colour_palette[4],
                'numerical_fisher_sp' : colour_palette[9],
                'numerical_fisher' : colour_palette[2],
                'heuristic' : colour_palette[1]}
