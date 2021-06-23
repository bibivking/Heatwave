# MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#  Author  : Mengyuan Mu
#  Version : 1.0 (22.06.2021)"
#  Email   : mu.mengyuan815@gmail.com
# WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

import matplotlib.pyplot as plt
from matplotlib.collections import EventCollection
import numpy as np


### Root distribution ###
PFT_number = [2, 5, 6, 9, 14, 2]	
rootbetas  = [0.962, 0.964, 0.943, 0.961, 0.961, 0.99]
PFT_names  = ['''evergreen 
broadleaf 
forest''',"shrub","grassland","cropland","bareland","modified"]

# Fixing random state for reproducibility
depths        = np.arange(0,460,10) 
root_fraction = np.zeros((len(PFT_number),len(depths))) 

for i,rootbeta in enumerate(rootbetas):
    for j,depth in enumerate(depths):
       root_fraction[i,j] = 1. - rootbeta**depth
print(root_fraction)

### Make plots ###
fig = plt.figure(figsize=(3.5,5))
fig.subplots_adjust(hspace=0.3)
fig.subplots_adjust(wspace=0.2)

plt.rcParams['text.usetex']     = False
plt.rcParams['font.family']     = "sans-serif"
plt.rcParams['font.serif']      = "Helvetica"
plt.rcParams['axes.linewidth']  = 1.5
plt.rcParams['axes.labelsize']  = 14
plt.rcParams['font.size']       = 14
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 14

almost_black = '#262626'
# change the tick colors also to the almost black
plt.rcParams['ytick.color']     = almost_black
plt.rcParams['xtick.color']     = almost_black

# change the text colors also to the almost black
plt.rcParams['text.color']      = almost_black

# Change the default axis colors from black to a slightly lighter black,
# and a little thinner (0.5 instead of 1)
plt.rcParams['axes.edgecolor']  = almost_black
plt.rcParams['axes.labelcolor'] = almost_black

# set the box type of sequence number
props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')

ax = fig.add_subplot(111)

# choose colormap
colors = [(0,0.5921569, 0.3607843),
          (0.7529412,0.8509804,0.5882353),
          (0.9882353,0.7647059,0.4666667),
          (0.9764706, 0.509804,0.2470588),
          (0.9333333,0.1568628,0.1254902),
          (0,0,1.)]
# (1,0.9882353,0.7568628)


for i,rootbeta in enumerate(rootbetas):
    ax.plot(root_fraction[i], depths, c=colors[i], lw=1.5, ls="-", label=PFT_names[i], alpha=1.)


ytickslocs     = [  0,    50,  100,   150,  200,   250,  300,   350,  400,  450]
yticklabels    = ["0.0", "0.5",  "1.0", "1.5",  "2.0", "2.5",  "3.0", "3.5",  "4.0", "4.5"]

xtickslocs     = [  0,  0.2,  0.4,  0.6,  0.8,    1.]
xticklabels    = ["0", "20", "40", "60", "80", "100"]

# plt.setp(ax.get_xticklabels(), visible=True)
ax.set(xticks=xtickslocs, xticklabels=xticklabels)
ax.set(yticks=ytickslocs, yticklabels=yticklabels)


ax.axis('tight')
ax.set_ylim(0.,460)
ax.set_ylabel("Depth (m)")
ax.set_xlim(0.,1.01)
ax.set_xlabel("Root Fraction (%)")
ax.invert_yaxis()

ax.text(0.02, 0.95, '(b)', transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
ax.legend(loc='lower left',prop={'size': 12}, frameon=False)

fig.savefig("./plots/FigS4_profile_root_distribution", bbox_inches='tight', pad_inches=0.1)

