import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib as mpl
from matplotlib import pylab
import numpy as np
from numpy import pi, sin, mean, nanmean, zeros, sqrt, multiply
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from numpy.fft import fft, fftfreq
from scipy import signal
import datetime
from matplotlib import dates as mdates
#############################
##########################################
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}
#matplotlib.rc('font', **font)
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
mpl.rcParams.update({'font.size': 20})

font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Lucida Grande']})
# for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
mpl.rcParams.update({'font.size': 20})
##########################################

#############################
'''
read and plot the w vs t
at z = -100 and plot spectra
'''

#############################
#############################
vlevel=-4500
suffix = 'velo_vlevel_' + str(vlevel)
suffix_full = 'velo_full_vlevel_' + str(vlevel)

nesm_surf_velo = np.loadtxt("nesm_avg_surf_velo_vlevel_" + str(vlevel) + ".asc")
# atlantis_surf_velo = np.loadtxt("atlantis_avg_surf_velo_vlevel_" + str(vlevel) + ".asc")
seamount_surf_velo = np.loadtxt("seamount_avg_surf_velo_vlevel_" + str(vlevel) + ".asc")

# nesm_surf_velo_full = np.loadtxt("nesm_avg_full_surf_velo_vlevel_" + str(vlevel) + ".asc")
# atlantis_surf_velo_full = np.loadtxt("atlantis_avg_full_surf_velo_vlevel_" + str(vlevel) + ".asc")
# seamount_surf_velo_full = np.loadtxt("seamount_avg_full_surf_velo_vlevel_" + str(vlevel) + ".asc")

base_date = datetime.datetime(2019, 1, 1)

# Create a date range
date_range = [base_date + datetime.timedelta(days=i) for i in range(len(nesm_surf_velo))]

#############################
linewidth=3
#############################
fig = plt.figure('avg_velo')
ax = fig.gca()  # Get current axes

ax.plot(date_range, nesm_surf_velo, linewidth=linewidth, label=r'\textrm{NESM}')
# ax.plot(date_range, atlantis_surf_velo, linewidth=linewidth, label=r'\textrm{ATLANTIS}')
ax.plot(date_range, seamount_surf_velo, linewidth=linewidth, label=r'\textrm{SEAMOUNT}')

# vlines_at_days = [90, 96, 140, 146, 263, 269, 378, 384]
vlines_at_days = [263, 269, 378, 384]

base_day = 0

for i in vlines_at_days:
        print("day in the current range = ", i - base_day)
        ax.axvline(x=date_range[i - base_day], color='b', linestyle='--', linewidth=0.5)

ax.set_xlabel(r"$t$",fontsize=18)
# ax.set_ylabel(r"\textrm{mean surface velocity}",fontsize=18)
ax.set_ylabel(r"$\left< {(u^2 + v^2)^{1/2}} \right>$",fontsize=18)
ax.set_xlim([date_range[41], date_range[-1]])


# Minor ticks every month.
fmt_month = mdates.MonthLocator()
# # Minor ticks every year.
fmt_year = mdates.YearLocator()

ax.xaxis.set_minor_locator(fmt_month)
# # '%b' to get the names of the month
ax.xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
ax.xaxis.set_major_locator(fmt_year)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

# # fontsize for month labels
ax.tick_params(labelsize=10, which='both')

# # create a second x-axis beneath the first x-axis to show the year in YYYY format
sec_xaxis = ax.secondary_xaxis(-0.1)
sec_xaxis.xaxis.set_major_locator(fmt_year)
sec_xaxis.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))


plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0)) 
plt.tight_layout()
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ratio = 0.5
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
ax.legend(loc=2,  prop={'size': 10})
fig.savefig('avg_' + suffix + '.pdf')
fig.savefig('avg_' + suffix + '.png')

#############################
#############################
# fig = plt.figure('avg_velo_full_hslice')
# ax = fig.gca()  # Get current axes

# ax.plot(nesm_surf_velo_full, linewidth=linewidth, label=r'\textrm{NESM}')
# ax.plot(atlantis_surf_velo_full, linewidth=linewidth, label=r'\textrm{ATLANTIS}')
# ax.plot(seamount_surf_velo_full, linewidth=linewidth, label=r'\textrm{SEAMOUNT}')
# # ax.axvline(x=65, color='k', linestyle='--', linewidth=0.5)
# ax.axvline(x=90, color='b', linestyle='--', linewidth=0.5)
# ax.axvline(x=96, color='b', linestyle='--', linewidth=0.5)

# # ax.axvline(x=120, color='b', linestyle='--', linewidth=0.5)
# # ax.axvline(x=126, color='b', linestyle='--', linewidth=0.5)

# ax.axvline(x=140, color='b', linestyle='--', linewidth=0.5)
# ax.axvline(x=146, color='b', linestyle='--', linewidth=0.5)

# ax.axvline(x=263, color='b', linestyle='--', linewidth=0.5)
# ax.axvline(x=269, color='b', linestyle='--', linewidth=0.5)

# ax.axvline(x=378, color='b', linestyle='--', linewidth=0.5)
# ax.axvline(x=384, color='b', linestyle='--', linewidth=0.5)

# ax.set_xlabel(r"$t$",fontsize=18)
# # ax.set_ylabel(r"\textrm{mean surface velocity}",fontsize=18)
# ax.set_ylabel(r"$\left< {(u^2 + v^2)^{1/2}} \right>$",fontsize=18)
# ax.set_xlim([41, 424])


# base = datetime.datetime(2019, 1, 1)
# # Minor ticks every month.
# fmt_month = mdates.MonthLocator()
# # Minor ticks every year.
# fmt_year = mdates.YearLocator()

# ax.xaxis.set_minor_locator(fmt_month)
# # '%b' to get the names of the month
# ax.xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
# ax.xaxis.set_major_locator(fmt_year)
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

# # # fontsize for month labels
# ax.tick_params(labelsize=10, which='both')

# # # create a second x-axis beneath the first x-axis to show the year in YYYY format
# # sec_xaxis = ax.secondary_xaxis(-0.1)
# # sec_xaxis.xaxis.set_major_locator(fmt_year)
# # sec_xaxis.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))


# plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0)) 
# plt.tight_layout()
# xleft, xright = ax.get_xlim()
# ybottom, ytop = ax.get_ylim()
# ratio = 0.5
# ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
# ax.legend(loc=2,  prop={'size': 10})
# fig.savefig('avg_' + suffix_full + '.pdf')
# fig.savefig('avg_' + suffix_full + '.png')

#############################

# find correlation coefficient
# corr_coef = np.corrcoef(nesm_surf_velo, atlantis_surf_velo)
# print("Correlation coefficient (nesm_surf_velo, atlantis_surf_velo):", corr_coef[0, 1])

corr_coef = np.corrcoef(nesm_surf_velo, seamount_surf_velo)
print("Correlation coefficient (nesm_surf_velo, seamount_surf_velo):", corr_coef[0, 1])

# corr_coef = np.corrcoef(atlantis_surf_velo, seamount_surf_velo)
# print("Correlation coefficient (atlantis_surf_velo, seamount_surf_velo):", corr_coef[0, 1])

#############################
surf_vlevel = -100
depth_vlevel = -4500
nesm_surf_velo = np.loadtxt("nesm_avg_surf_velo_vlevel_" + str(surf_vlevel) + ".asc")
seamount_surf_velo = np.loadtxt("seamount_avg_surf_velo_vlevel_" + str(surf_vlevel) + ".asc")
nesm_depth_velo = np.loadtxt("nesm_avg_surf_velo_vlevel_" + str(depth_vlevel) + ".asc")
seamount_depth_velo = np.loadtxt("seamount_avg_surf_velo_vlevel_" + str(depth_vlevel) + ".asc")

corr_coef = np.corrcoef(nesm_surf_velo, nesm_depth_velo)
print("Correlation coefficient (nesm_surf_velo, nesm_depth_velo):", corr_coef[0, 1])

corr_coef = np.corrcoef(seamount_surf_velo, seamount_depth_velo)
print("Correlation coefficient (seamount_surf_velo, seamount_depth_velo):", corr_coef[0, 1])

corr_coef = np.corrcoef(nesm_depth_velo, seamount_depth_velo)
print("Correlation coefficient (nesm_depth_velo, seamount_depth_velo):", corr_coef[0, 1])
