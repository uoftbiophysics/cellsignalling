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
KON = 1E5
KP = 1.
T = 1E2
KF = None
ALPHA = 0.2
C1 = 1E-6 # concentration nano to micro molar
C2 = 1E-6
KOFF = 1E-2 # roughly 10^-2
KOFF2 = 1E0
N = 1E2

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
