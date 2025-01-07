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
suffix = '_t'
# suffix_full = 'velo_full_vlevel_' + str(vlevel)
vlevel = -4500

nesm_kmke = np.loadtxt("nesm_2019_2020_KmKe_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")
atlantis_kmke = np.loadtxt("atlantis_2019_2020_KmKe_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")
seamount_kmke = np.loadtxt("seamount_2019_2020_KmKe_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")

nesm_hrs = np.loadtxt("nesm_2019_2020_HRS_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")
atlantis_hrs = np.loadtxt("atlantis_2019_2020_HRS_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")
seamount_hrs = np.loadtxt("seamount_2019_2020_HRS_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")

nesm_vrs = np.loadtxt("nesm_2019_2020_VRS_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")
atlantis_vrs = np.loadtxt("atlantis_2019_2020_VRS_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")
seamount_vrs = np.loadtxt("seamount_2019_2020_VRS_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")

nesm_vebf = np.loadtxt("nesm_2019_2020_VEBF_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")
atlantis_vebf = np.loadtxt("atlantis_2019_2020_VEBF_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")
seamount_vebf = np.loadtxt("seamount_2019_2020_VEBF_nt_" + str(indxRange) + "_vlevel_" + str(vlevel) + ".asc")

#############################
linewidth=1
linewidth_hrs = 1
linewidth_vrs = 1
nesm_kmke_linecolor = 'tab:blue'
nesm_vebf_linecolor = 'tab:blue'

# atlantis_kmke_linecolor = 'tab:orange'
# atlantis_vebf_linecolor = 'tab:orange'

seamount_kmke_linecolor = 'tab:orange'
seamount_vebf_linecolor = 'tab:orange'

hrs_linecolor = 'k'
vrs_linecolor = 'k'

kmke_marker = 'none'
vebf_marker = 'none'

markersize=5
markerskip=5

kmke_min = -1.2e-7
kmke_max = 1.2e-7
vebf_min = -1.2e-7;
vebf_max = 1.2e-7
#############################
Nt = len(nesm_kmke)
n_files_per_day = 8
Nt_avg = int(Nt/n_files_per_day)
print("Nt = ", Nt, "n_files_per_day = ", n_files_per_day, \
        "Nt_avg = int(Nt/n_files_per_day) =", Nt_avg)

nesm_kmke_daily = np.zeros(Nt_avg)
atlantis_kmke_daily = np.zeros(Nt_avg)
seamount_kmke_daily = np.zeros(Nt_avg)

nesm_hrs_daily = np.zeros(Nt_avg)
atlantis_hrs_daily = np.zeros(Nt_avg)
seamount_hrs_daily = np.zeros(Nt_avg)

nesm_vrs_daily = np.zeros(Nt_avg)
atlantis_vrs_daily = np.zeros(Nt_avg)
seamount_vrs_daily = np.zeros(Nt_avg)

nesm_vebf_daily = np.zeros(Nt_avg)
atlantis_vebf_daily = np.zeros(Nt_avg)
seamount_vebf_daily = np.zeros(Nt_avg)

# time avg
for ll in range(0, Nt_avg):
        nesm_kmke_daily[ll] =  (nanmean(nesm_kmke[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        atlantis_kmke_daily[ll] =  (nanmean(atlantis_kmke[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        seamount_kmke_daily[ll] =  (nanmean(seamount_kmke[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))

        nesm_hrs_daily[ll] =  (nanmean(nesm_hrs[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        atlantis_hrs_daily[ll] =  (nanmean(atlantis_hrs[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        seamount_hrs_daily[ll] =  (nanmean(seamount_hrs[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))

        nesm_vrs_daily[ll] =  (nanmean(nesm_vrs[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        atlantis_vrs_daily[ll] =  (nanmean(atlantis_vrs[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        seamount_vrs_daily[ll] =  (nanmean(seamount_vrs[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))

        nesm_vebf_daily[ll] =  (nanmean(nesm_vebf[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        atlantis_vebf_daily[ll] =  (nanmean(atlantis_vebf[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
        seamount_vebf_daily[ll] =  (nanmean(seamount_vebf[ll*n_files_per_day:(ll+1)*(n_files_per_day)]))

base_date = datetime.datetime(2019, 3, 1)

# Create a date range
date_range = [base_date + datetime.timedelta(days=i) for i in range(len(nesm_kmke_daily))]

# print(date_range)
#############################
fig = plt.figure('kmke')
ax = fig.gca()  # Get current axes

ax.plot(date_range, nesm_kmke_daily, \
        color=nesm_kmke_linecolor, linewidth=linewidth, \
        marker=kmke_marker, markersize=markersize, markevery=markerskip, \
        label=r'\textrm{NESM}')
# ax.plot(date_range, atlantis_kmke_daily, \
#         color=atlantis_kmke_linecolor, linewidth=linewidth, \
#         marker=kmke_marker, markersize=markersize, markevery=markerskip, \
#         label=r'\textrm{ATLANTIS}')
ax.plot(date_range, seamount_kmke_daily, \
        color=seamount_kmke_linecolor, linewidth=linewidth, \
        marker=kmke_marker, markersize=markersize, markevery=markerskip, \
        label=r'\textrm{SEAMOUNT}')

# ax.plot(date_range, nesm_hrs_daily, \
#         color=hrs_linecolor, linewidth=linewidth_hrs, label=r'\textrm{NESM} HRS')
# ax.plot(date_range, atlantis_hrs_daily, \
#         color=hrs_linecolor, linewidth=linewidth_hrs, label=r'\textrm{ATLANTIS} HRS')
# ax.plot(date_range, seamount_hrs_daily, \
#         color=hrs_linecolor, linewidth=linewidth_hrs, label=r'\textrm{SEAMOUNT} HRS')

# ax.plot(date_range, nesm_vrs_daily, \
#         color=vrs_linecolor, linewidth=linewidth_vrs, label=r'\textrm{NESM} VRS')
# ax.plot(date_range, atlantis_vrs_daily, \
#         color=vrs_linecolor, linewidth=linewidth_vrs, label=r'\textrm{ATLANTIS} VRS')
# ax.plot(date_range, seamount_vrs_daily, \
#         color=vrs_linecolor, linewidth=linewidth_vrs, label=r'\textrm{SEAMOUNT} VRS')

# vlines_at_days = [90, 96, 140, 146, 263, 269, 378, 384]
vlines_at_days = [263, 269, 378, 384]
base_day = 60

for i in vlines_at_days:
        print("day in the current range = ", i - base_day)
        ax.axvline(x=date_range[i - base_day], color='b', linestyle='--', linewidth=0.5)

ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
ax.set_xlabel(r"$t$",fontsize=18)
# ax.set_ylabel(r"\textrm{mean surface velocity}",fontsize=18)
ax.set_ylabel(r"$\left< KmKe \right>$",fontsize=18)
ax.set_xlim([date_range[0], date_range[-1]])
ax.set_ylim([kmke_min, kmke_max])


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
fig.savefig('avg_kmke' + suffix + '.pdf')
fig.savefig('avg_kmke' + suffix + '.png')

#############################

fig = plt.figure('vebf')
ax = fig.gca()  # Get current axes

ax.plot(date_range, nesm_vebf_daily, \
        color=nesm_vebf_linecolor, linewidth=linewidth, \
        marker=vebf_marker, markersize=markersize, markevery=markerskip, \
        label=r'\textrm{NESM}')
# ax.plot(date_range, atlantis_vebf_daily, \
#         color=atlantis_vebf_linecolor, linewidth=linewidth, \
#         marker=vebf_marker, markersize=markersize, markevery=markerskip, \
#         label=r'\textrm{ATLANTIS}')
ax.plot(date_range, seamount_vebf_daily, \
        color=seamount_vebf_linecolor, linewidth=linewidth, \
        marker=vebf_marker, markersize=markersize, markevery=markerskip, \
        label=r'\textrm{SEAMOUNT}')

# vlines_at_days = [90, 96, 140, 146, 263, 269, 378, 384]
vlines_at_days = [263, 269, 378, 384]

base_day = 60

for i in vlines_at_days:
        print("day in the current range = ", i - base_day)
        ax.axvline(x=date_range[i - base_day], color='b', linestyle='--', linewidth=0.5)

ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)

ax.set_xlabel(r"$t$",fontsize=18)
# ax.set_ylabel(r"\textrm{mean surface velocity}",fontsize=18)
ax.set_ylabel(r"$\left< VEBF \right>$",fontsize=18)
ax.set_xlim([date_range[0], date_range[-1]])
ax.set_ylim([vebf_min, vebf_max])



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
fig.savefig('avg_vebf' + suffix + '.pdf')
fig.savefig('avg_vebf' + suffix + '.png')

#############################
