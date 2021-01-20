import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
print(sns.__version__)

fname = "./txt/scatter_WTD_Dr_CTL_PFT-tree_2010-11.txt"
x = pd.read_csv(fname, names=["wt_depth","delta_drain"], sep=" ", skipinitialspace=True)
sns.displot(x, x="wt_depth", y="delta_drain")#, xlim=[0,10], ylim=[0,15])
# plt.xlim(0, 15)
# plt.ylim(10, 100)
plt.savefig("scatter_density_Dr_2010-11.png")

#
# fig = plt.figure(figsize=[15,10])
#
# ax  = fig.add_subplot(111)
#
# fig.subplots_adjust(hspace=0.1)
# fig.subplots_adjust(wspace=0.05)
# plt.rcParams['text.usetex'] = False
# plt.rcParams['font.family'] = "sans-serif"
# plt.rcParams['font.sans-serif'] = "Helvetica"
# plt.rcParams['axes.labelsize'] = 14
# plt.rcParams['font.size'] = 14
# plt.rcParams['legend.fontsize'] = 12
# plt.rcParams['xtick.labelsize'] = 14
# plt.rcParams['ytick.labelsize'] = 14
#
# almost_black = '#262626'
# # change the tick colors also to the almost black
# plt.rcParams['ytick.color'] = almost_black
# plt.rcParams['xtick.color'] = almost_black
#
# # change the text colors also to the almost black
# plt.rcParams['text.color'] = almost_black
#
# # Change the default axis colors from black to a slightly lighter black,
# # and a little thinner (0.5 instead of 1)
# plt.rcParams['axes.edgecolor'] = almost_black
# plt.rcParams['axes.labelcolor'] = almost_black
# ax.plot(x, fw1["cable"],   c="orange", lw=1.0, ls="-", label="Default_fw-std")
# ax.plot(x, fw2["cable"],   c="blue", lw=1.0, ls="-", label="Best_fw-std")
# ax.plot(x, fw3["cable"],   c="forestgreen", lw=1.0, ls="-", label="Best_fw-hie")
#
# cleaner_dates = ["2013","2014","2015","2016","2017","2018","2019"]
# xtickslocs    = [367,732,1097,1462,1828,2193,2558]
#
# plt.setp(ax.get_xticklabels(), visible=True)
# ax.set(xticks=xtickslocs, xticklabels=cleaner_dates)
# ax.set_ylabel("β")
# ax.axis('tight')
# ax.set_ylim(0.,1.1)
# ax.set_xlim(367,2739)
# ax.legend()
# fig.savefig("../plots/EucFACE_fwsoil_comp_%s" % ring, bbox_inches='tight', pad_inches=0.1)



# fig = sns.histplot(x, x="wt_depth", y="Qrec")
# fig = sns.distplot(x, x="wt_depth", y="Qrec")
# fig = sns.jointplot(data=x, x="wt_depth", y="Qrec", xlim=[0,10], ylim=[0,15], alpha=0.2, kind="scatter")
# fig.savefig("scatter_density_drainage.png")

# df = pd.read_csv(fname, names=['wt_depth', 'Qrec'], sep=' ', skipinitialspace=True)
# sns.histplot(df, x="wt_depth", y="Qrec", xlim=[0,10], ylim=[0,15], alpha=0.2, )
# sns.plt.show()
