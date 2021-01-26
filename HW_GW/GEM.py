import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pygam import LinearGAM
# from matplotlib import cm

#from pygam.utils import generate_X_grid
fname = "./txt/scatter_WTD_deltaT_CTL_PFT-tree_2000-19.txt"
df = pd.read_csv(fname, names=["WTD (m)", "ΔT (°C)"], sep=' ', skipinitialspace=True)

width  = 6
height = 6
fig    = plt.figure(figsize=(width, height))
fig.subplots_adjust(hspace=0.1)
fig.subplots_adjust(wspace=0.05)
plt.rcParams['text.usetex']     = False
plt.rcParams['font.family']     = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize']  = 18
plt.rcParams['font.size']       = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

# props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
# colors = cm.Oranges()
# , c=colors,
print("I am OK 1 ")
ax = fig.add_subplot(111)

sns.histplot(
    df, x="WTD (m)", y="ΔT (°C)",
    bins=(100,100), discrete=(False, False),
    cbar=False, cbar_kws=dict(shrink=.75), ax = ax)
# sns.displot(df, x="WTD", y="deltaT",ax = ax)
x = df['WTD (m)'].values
y = df['ΔT (°C)'].values
xx = x.reshape(x.shape[0], 1)
yy = y.reshape(y.shape[0], 1)

print("I am OK 2")
# reshape for gam
gam    = LinearGAM(n_splines=20).gridsearch(xx, yy)
x_pred = np.linspace(min(x), max(x), num=100)
y_pred = gam.predict(x_pred)
y_int  = gam.confidence_intervals(x_pred, width=.95)
ax.plot(x_pred, y_pred, color="red", ls='-', lw=2.0, zorder=10)
ax.fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
                facecolor='red', zorder=10)
# ax.text(0.03, 0.95, '(f)', transform=ax.transAxes, fontsize=18, verticalalignment='top', bbox=props)

print("I am OK 3")
# plt.setp(ax.get_xticklabels(), visible=False)
# ax.set(xticks=xtickslocs, xticklabels=cleaner_dates) ####
ax.yaxis.tick_left()
ax.yaxis.set_label_position("left")
# ax.set_ylabel("Evaporative Fraction (-)")
# ax.axis('tight')
#ax.set_ylim(0.,120.)
#ax.set_xlim(367,2739)#,1098)
# ax.set_xlim(day_start,2739)

# ax.set_ylabel("Top soil moisture  (m$3$ m$-3$)")
# ax.axis('tight')

ax.set_xlim(0,10)
ax.set_ylim(-4,8)
fig.savefig("scatter_density_deltaT_fitting_2000-19.pdf")
