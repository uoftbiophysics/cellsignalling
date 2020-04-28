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

# globals for heatmaps and equations (to prevent circular imports)
KON = 1.
KP = 100.
T = 100.
KF = 1.0
ALPHA = 0.2
C1 = 10.0
C2 = 1.0
KOFF = 1.0
KOFF2 = 10.0

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
