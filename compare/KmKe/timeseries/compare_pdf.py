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
from scipy.stats import gaussian_kde

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
suffix = '_pdf'
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
# Create a Gaussian KDE
nesm_kmke_kde = gaussian_kde(nesm_kmke_daily)
atlantis_kmke_kde = gaussian_kde(atlantis_kmke_daily)
seamount_kmke_kde = gaussian_kde(seamount_kmke_daily)

nesm_vebf_kde = gaussian_kde(nesm_vebf_daily)
atlantis_vebf_kde = gaussian_kde(atlantis_vebf_daily)
seamount_vebf_kde = gaussian_kde(seamount_vebf_daily)


# Create an array of values over which to evaluate the KDE
kmke_x = np.linspace(kmke_min, kmke_max, 1000)
vebf_x = np.linspace(vebf_min, vebf_max, 1000)

# Evaluate the PDF
nesm_kmke_pdf = nesm_kmke_kde(kmke_x)
# nesm_kmke_pdf = nesm_kmke_pdf/np.trapz(nesm_kmke_pdf, kmke_x)

atlantis_kmke_pdf = atlantis_kmke_kde(kmke_x)
# atlantis_kmke_pdf = atlantis_kmke_pdf/np.trapz(atlantis_kmke_pdf, kmke_x)

seamount_kmke_pdf = seamount_kmke_kde(kmke_x)
# seamount_kmke_pdf = seamount_kmke_pdf/np.trapz(seamount_kmke_pdf, kmke_x)

nesm_vebf_pdf = nesm_vebf_kde(vebf_x)
# nesm_vebf_pdf = nesm_vebf_pdf/np.trapz(nesm_vebf_pdf, vebf_x)

atlantis_vebf_pdf = atlantis_vebf_kde(vebf_x)
# atlantis_vebf_pdf = atlantis_vebf_pdf/np.trapz(atlantis_vebf_pdf, vebf_x)

seamount_vebf_pdf = seamount_vebf_kde(vebf_x)
# seamount_vebf_pdf = seamount_vebf_pdf/np.trapz(seamount_vebf_pdf, vebf_x)


#############################
#############################
fig = plt.figure('kmke')
ax = fig.gca()  # Get current axes

ax.plot(kmke_x, nesm_kmke_pdf, label=r'\textrm{NESM}')
# ax.plot(kmke_x, atlantis_kmke_pdf, label=r'\textrm{ATLANTIS}')
ax.plot(kmke_x, seamount_kmke_pdf, label=r'\textrm{SEAMOUNT}')

ax.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
ax.set_xlabel(r"$\left< KmKe \right>$",fontsize=18)
# ax.set_ylabel(r"\textrm{mean surface velocity}",fontsize=18)
ax.set_ylabel(r"$PDF$",fontsize=18)

xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ratio = 0.5
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
ax.legend(loc=2,  prop={'size': 10})
fig.savefig('avg_kmke' + suffix + '.pdf')
fig.savefig('avg_kmke' + suffix + '.png')

#############################
#############################
fig = plt.figure('vebf')
ax = fig.gca()  # Get current axes

ax.plot(vebf_x, nesm_vebf_pdf, label=r'\textrm{NESM}')
# ax.plot(vebf_x, atlantis_vebf_pdf, label=r'\textrm{ATLANTIS}')
ax.plot(vebf_x, seamount_vebf_pdf, label=r'\textrm{SEAMOUNT}')

ax.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
ax.set_xlabel(r"$\left< VEBF \right>$",fontsize=18)
# ax.set_ylabel(r"\textrm{mean surface velocity}",fontsize=18)
ax.set_ylabel(r"$PDF$",fontsize=18)

xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ratio = 0.5
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
ax.legend(loc=2,  prop={'size': 10})
fig.savefig('avg_vebf' + suffix + '.pdf')
fig.savefig('avg_vebf' + suffix + '.png')

#############################
