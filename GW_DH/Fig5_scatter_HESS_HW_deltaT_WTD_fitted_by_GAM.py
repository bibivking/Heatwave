# MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#  Author  : Mengyuan Mu
#  Version : 1.0 (22.03.2021)"
#  Email   : mu.mengyuan815@gmail.com
# WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pygam import LinearGAM
# from matplotlib import cm

#from pygam.utils import generate_X_grid
fname = "./txt/scatter_WTD_deltaT_CTL_PFT-tree_2000-19.txt"
df = pd.read_csv(fname, names=["WTD (m)", "ΔT (°C)"], sep=' ', skipinitialspace=True)

# props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
# colors = cm.Oranges()

print("I am OK 1 ")

'''
    width  = 7
    height = 7
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

    # histplot
    ax = fig.add_subplot(111)
    # sns.distplot(df,x="WTD (m)", y="ΔT (°C)",cbar=False, cbar_kws=dict(shrink=.75))
    sns.histplot(
        df, x="WTD (m)", y="ΔT (°C)",
        bins=(100,100), discrete=(False, False),
        cbar=False, cbar_kws=dict(shrink=.75), ax = ax)

    x = df['WTD (m)'].values
    y = df['ΔT (°C)'].values
    xx = x.reshape(x.shape[0], 1)
    yy = y.reshape(y.shape[0], 1)

    # reshape for gam
    # gam    = LinearGAM(n_splines=4).gridsearch(xx, yy) # n_splines=22
    # x_pred = np.linspace(min(x), max(x), num=100)
    # y_pred = gam.predict(x_pred)
    # y_int  = gam.confidence_intervals(x_pred, width=.95)
    # ax.plot(x_pred, y_pred, color="red", ls='-', lw=2.0, zorder=10)
    # ax.fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
    #                 facecolor='red', zorder=10)
    # ax.text(0.03, 0.95, '(f)', transform=ax.transAxes, fontsize=18, verticalalignment='top', bbox=props)

    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")
    bwith = 2
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    ax.set_xlim(0,10)
    ax.set_ylim(-4,8)
    fig.savefig("./plots/Fig5_scatter_density_deltaT_fitting_2000-19.pdf")
'''



# displot
# width  = 7
# height = 7
# fig    = plt.figure(figsize=(width, height))
# fig.subplots_adjust(hspace=0.1)
# fig.subplots_adjust(wspace=0.05)


g = sns.displot(df,x="WTD (m)", y="ΔT (°C)",cbar=False, cbar_kws=dict(shrink=.75),
                height=3, aspect=1.,facet_kws=dict(despine=False))
g.set(xlim=(0, 10), ylim=(-4, 8))
# g.set(font=dict(family="sans-serif",size=18),
#       axes=dict(labelsize=18),legend=dict(fontsize=18),
#       xtick=dict(labelsize=18),ytick=dict(labelsize=18))
#text=dict(usetex=False),
g.savefig("./plots/Fig5_scatter_density_deltaT_fitting_2000-19.pdf")
# scatterplot
# f, ax = plt.subplots(figsize=(7, 7))
# sns.despine(f, left=True, bottom=True)
# sns.scatterplot(data=df, x="WTD (m)", y="ΔT (°C)",ax=ax)
# f.savefig("./plots/Fig5_scatter_density_deltaT_fitting_2000-19.pdf")
#
