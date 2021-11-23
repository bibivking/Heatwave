#!/usr/bin/env python
"""
Produce the Australia map with the label of EucFACE site - Figure 1
"""

__author__ = "Mengyuan Mu & Martin De Kauwe"
__email__  = "mu.mengyuan815@gmail.com"

import os
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from shapely.geometry import LineString, MultiLineString
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import sys
import matplotlib.ticker as mticker
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader

fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.size'] = "14"
plt.rcParams['font.sans-serif'] = "Helvetica"

ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.COASTLINE, lw=0.5)

ax.set_extent([-150, 60, -25, 60])

from cartopy.feature import NaturalEarthFeature

coast = NaturalEarthFeature(category='physical', scale='10m',
                            facecolor='none', name='coastline')
feature = ax.add_feature(coast, edgecolor='black', lw=0.5)

import cartopy.feature as cfeature

# add and color high-resolution land
LAND_highres = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='white',
                                            facecolor='#fff9e6',
                                            linewidth=.1
                                           )
ax.add_feature(LAND_highres, zorder=0,
               edgecolor='#fff9e6', facecolor='#fff9e6')

ax.set_xlim(101, 169)
ax.set_ylim(-49, -1)

ax.plot(151.2093, -33.8688, 'ko', markersize=4, transform=ccrs.PlateCarree())
ax.text(155.2, -32, 'Sydney', horizontalalignment='right',
        transform=ccrs.PlateCarree(), fontsize=14)

ax.text(141, -26, 'Australia', horizontalalignment='right',
        transform=ccrs.PlateCarree(), fontsize=20)

ax.text(110, -5, '(a)', horizontalalignment='right',
        transform=ccrs.PlateCarree(), fontsize=16)

# Add EucFACE as a star
ax.plot(150.738, -33.62, '*', markersize=20, c= "red", markeredgecolor = 'red', alpha = 0.7,
        transform=ccrs.PlateCarree())

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='black', alpha=0.5,
                  linestyle='--')

gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = True
gl.ylines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

gl.xlocator = mticker.FixedLocator([100, 110, 120, 130, 140, 150, 160, 170])
gl.ylocator = mticker.FixedLocator([0, -10, -20, -30, -40, -50])

fdir = "./plots"
fig.savefig(os.path.join(fdir, "aus_map.png"), dpi=300, bbox_inches='tight',
            pad_inches=0.1)
© 2021 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About

