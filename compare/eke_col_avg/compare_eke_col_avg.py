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
# vlevel=-4500
indxRange = '469_3382'
suffix = 'eke_col_avg'
# suffix_full = 'velo_full_vlevel_' + str(vlevel)

nesm_surf_velo = np.loadtxt("nesm_2019_2020_eke_t_col_avg_nt_" + str(indxRange) + ".asc")
# atlantis_surf_velo = np.loadtxt("atlantis_2019_2020_eke_t_col_avg_nt_" + str(indxRange) + ".asc")
seamount_surf_velo = np.loadtxt("seamount_2019_2020_eke_t_col_avg_nt_" + str(indxRange) + ".asc")



#############################
linewidth=3

Nt = len(nesm_surf_velo)
n_files_per_day = 8
Nt_avg = int(Nt/n_files_per_day)
print("Nt = ", Nt, "n_files_per_day = ", n_files_per_day, \
        "Nt_avg = int(Nt/n_files_per_day) =", Nt_avg)

nesm_eke = np.zeros(Nt_avg)
atlantis_eke = np.zeros(Nt_avg)
seamount_eke = np.zeros(Nt_avg)
# time avg
for ll in range(0, Nt_avg):
        nesm_eke[ll] =  (nanmean(nesm_surf_velo[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        # atlantis_eke[ll] =  (nanmean(atlantis_surf_velo[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        seamount_eke[ll] =  (nanmean(seamount_surf_velo[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))

base_date = datetime.datetime(2019, 3, 1)

# Create a date range
date_range = [base_date + datetime.timedelta(days=i) for i in range(len(nesm_eke))]

# print(date_range)
#############################
fig = plt.figure('eke')
ax = fig.gca()  # Get current axes

ax.plot(date_range, nesm_eke, linewidth=linewidth, label=r'\textrm{NESM}')
# ax.plot(date_range, atlantis_eke, linewidth=linewidth, label=r'\textrm{ATLANTIS}')
ax.plot(date_range, seamount_eke, linewidth=linewidth, label=r'\textrm{SEAMOUNT}')

# vlines_at_days = [90, 96, 140, 146, 263, 269, 378, 384]
vlines_at_days = [263, 269, 378, 384]

base_day = 60

for i in vlines_at_days:
        print("day in the current range = ", i - base_day)
        ax.axvline(x=date_range[i - base_day], color='b', linestyle='--', linewidth=0.5)


ax.set_xlabel(r"$t$",fontsize=18)
# ax.set_ylabel(r"\textrm{mean surface velocity}",fontsize=18)
ax.set_ylabel(r"$\left< EKE \right>$",fontsize=18)
ax.set_xlim([date_range[0], date_range[-1]])



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

# find correlation coefficient
# corr_coef = np.corrcoef(nesm_surf_velo, atlantis_surf_velo)
# print("Correlation coefficient (nesm_surf_velo, atlantis_surf_velo):", corr_coef[0, 1])

# corr_coef = np.corrcoef(nesm_surf_velo, seamount_surf_velo)
# print("Correlation coefficient (nesm_surf_velo, seamount_surf_velo):", corr_coef[0, 1])

# corr_coef = np.corrcoef(atlantis_surf_velo, seamount_surf_velo)
# print("Correlation coefficient (atlantis_surf_velo, seamount_surf_velo):", corr_coef[0, 1])

