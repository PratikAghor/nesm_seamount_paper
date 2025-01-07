import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib as mpl
from matplotlib import pylab
import numpy as np
from numpy import pi, sin, mean, nanmean, zeros, sqrt, multiply,shape
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from numpy.fft import fft, fftfreq
from scipy import signal
#############################
"""
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
"""
#############################
'''
read and plot the w vs t
at z = -1000 and plot spectra
'''

#############################
# indxRange='469_3382' # entire year
# indxRange = '712_768'; # Apr 01 - Apr 07, 2019
# indxRange = '1040_1096'; # May 12 - May 19, 2019
# indxRange='2104_2160' # Sept 22 - Sept 29, 2020
indxRange = '3016_3072'; # Jan 14 - Jan 20, 2020
# indxRange='3080_3144' # Jan 2020
lat_arr = np.array([38.3049, 38.5027, 38.7095])

lon_arr = np.array([-63.6058, -63.3725, -63.0808, -62.9058])


if(indxRange=='469_3382'):
        max_KmKe_ylim = 5e-7
        min_KmKe_ylim = -5e-7
else:
        max_KmKe_ylim = 1e-6
        min_KmKe_ylim = -1e-6


markerskip=5
markersize=5
#############################
# box_arr = np.array(["N", "NE", "W", "SW", "S", "SE", "NW", "E"])
box_arr =  np.array(["big_N", "big_S"])
for i in range (0,len(box_arr)):
        box = box_arr[i]
        
        suffix = box + '_nt_' + str(indxRange)
        nesm_hrs = np.loadtxt('hrs/nesm_2019_2020_hrs_z_'+ box + '_nt_'+ str(indxRange) + '.asc')
        # atlantis_hrs = np.loadtxt('hrs/atlantis_2019_2020_hrs_z_'+ box +'_nt_'+ str(indxRange) + '.asc')
        seamount_hrs = np.loadtxt('hrs/seamount_2019_2020_hrs_z_' + box + '_nt_'+ str(indxRange) + '.asc')

        nesm_vrs = np.loadtxt('vrs/nesm_2019_2020_vrs_z_'+ box +'_nt_'+ str(indxRange) + '.asc')
        # atlantis_vrs = np.loadtxt('vrs/atlantis_2019_2020_vrs_z_'+ box +'_nt_'+ str(indxRange) + '.asc')
        seamount_vrs = np.loadtxt('vrs/seamount_2019_2020_vrs_z_' + box + '_nt_'+ str(indxRange) + '.asc')

        nesm_kmke = np.loadtxt('kmke/nesm_2019_2020_kmke_z_'+ box +'_nt_'+ str(indxRange) + '.asc')
        # atlantis_kmke = np.loadtxt('kmke/atlantis_2019_2020_kmke_z_'+ box +'_nt_'+ str(indxRange) + '.asc')
        seamount_kmke = np.loadtxt('kmke/seamount_2019_2020_kmke_z_' + box + '_nt_'+ str(indxRange) + '.asc')

        nesm_vebf = np.loadtxt('vebf/nesm_2019_2020_vebf_z_' + box + '_nt_'+ str(indxRange) + '.asc')
        # atlantis_vebf = np.loadtxt('vebf/atlantis_2019_2020_vebf_z_' + box + '_nt_'+ str(indxRange) + '.asc')
        seamount_vebf = np.loadtxt('vebf/seamount_2019_2020_vebf_z_' + box + '_nt_'+ str(indxRange) + '.asc')
        '''
        nesm_kmke = np.zeros(np.shape(nesm_hrs))
        atlantis_kmke = np.zeros(np.shape(atlantis_hrs))
        seamount_kmke = np.zeros(np.shape(seamount_hrs))

        nesm_kmke[:, 0] = nesm_hrs[:, 0]  # z values
        nesm_kmke[:, 1] = nesm_hrs[:, 1] + nesm_vrs[:, 1] # kmke = hrs + vrs

        atlantis_kmke[:, 0] = atlantis_hrs[:, 0]  # z values
        atlantis_kmke[:, 1] = atlantis_hrs[:, 1] + atlantis_vrs[:, 1] # kmke = hrs + vrs

        seamount_kmke[:, 0] = seamount_hrs[:, 0]  # z values
        seamount_kmke[:, 1] = seamount_hrs[:, 1] + seamount_vrs[:, 1] # kmke = hrs + vrs
        ''' 
        #############################

        #############################
        fig = plt.figure('nesm_kmke_decomp_')
        ax = fig.gca()  # Get current axes
        #-----------------------------------------
        # ax.plot(seamount_hrs[:, 0], seamount_hrs[:, 1], \
        #        linewidth=1, linestyle='--', marker='>', markersize=markersize, \
        #        markevery=markerskip, color='k', label=r'\textrm{HRS}')
        #-----------------------------------------
        #ax.plot(seamount_vrs[:, 0], seamount_vrs[:, 1], \
        #        linewidth=1, linestyle='--', marker='s', markersize=markersize, \
        #        markevery=markerskip, color='k', label=r'\textrm{VRS}')
        #-----------------------------------------
        #ax.plot(seamount_kmke[:, 0], seamount_kmke[:, 1], \
        #        linewidth=2, linestyle='-', marker='o', markersize=markersize, \
        #        markevery=markerskip, color='r', label=r'\textrm{KmKe}')        
        # ax.set_xlabel(r"$z$",fontsize=18)
        # ax.set_ylabel(r"$KmKe$",fontsize=18)


        # without latex labels
        #-----------------------------------------
        ax.plot(nesm_hrs[:, 0], nesm_hrs[:, 1], \
                linewidth=1, linestyle='--', marker='>', markersize=markersize, \
                markevery=markerskip, color='k', label='HRS')
        #-----------------------------------------
        ax.plot(nesm_vrs[:, 0], nesm_vrs[:, 1], \
                linewidth=1, linestyle='--', marker='s', markersize=markersize, \
                markevery=markerskip, color='k', label='VRS')
        #-----------------------------------------
        ax.plot(nesm_kmke[:, 0], nesm_kmke[:, 1], \
                linewidth=2, linestyle='-', marker='o', markersize=markersize, \
                markevery=markerskip, color='r', label='KmKe')  
        #-----------------------------------------
        ax.plot(nesm_vebf[:, 0], nesm_vebf[:, 1], \
                linewidth=2, linestyle='-', marker='x', markersize=markersize, \
                markevery=markerskip, color='b', label='VEBF')
        #-----------------------------------------        
        ax.set_xlabel("z",fontsize=18)
        ax.set_ylabel("KmKe",fontsize=18)
        ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)

        ax.set_xlim([0, -5000])
        ax.set_ylim([min_KmKe_ylim, max_KmKe_ylim])

        xticks= np.array([-5000, -4000, -3000, -2000, -1000,])
        # xticklabels = np.array([r'$5000$', r'$4000$', r'$3000$', r'$2000$', r'$1000$'])
        xticklabels = np.array(['5000', '4000', '3000', '2000', '1000'])
        ax.set_xticks(xticks) 
        ax.set_xticklabels(xticklabels)

        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0)) 
        plt.tight_layout()
        xleft, xright = ax.get_xlim()
        ybottom, ytop = ax.get_ylim()
        ratio = 0.5
        ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
        ax.legend(loc='best',  prop={'size': 10})
        fig.savefig(str(indxRange) +'/nesm_avg_kmke_decomp_' + suffix + '.pdf')
        # fig.savefig(str(indxRange) + '/nesm_avg_kmke_decomp_' + suffix + '.png')
        plt.close() # clear figures to avoid overwriting in loop
        #############################
        #############################
        fig = plt.figure('seamount_kmke_decomp_')
        ax = fig.gca()  # Get current axes
        #-----------------------------------------
        # ax.plot(seamount_hrs[:, 0], seamount_hrs[:, 1], \
        #        linewidth=1, linestyle='--', marker='>', markersize=markersize, \
        #        markevery=markerskip, color='k', label=r'\textrm{HRS}')
        #-----------------------------------------
        #ax.plot(seamount_vrs[:, 0], seamount_vrs[:, 1], \
        #        linewidth=1, linestyle='--', marker='s', markersize=markersize, \
        #        markevery=markerskip, color='k', label=r'\textrm{VRS}')
        #-----------------------------------------
        #ax.plot(seamount_kmke[:, 0], seamount_kmke[:, 1], \
        #        linewidth=2, linestyle='-', marker='o', markersize=markersize, \
        #        markevery=markerskip, color='r', label=r'\textrm{KmKe}')        
        # ax.set_xlabel(r"$z$",fontsize=18)
        # ax.set_ylabel(r"$KmKe$",fontsize=18)


        # without latex labels
        #-----------------------------------------
        ax.plot(seamount_hrs[:, 0], seamount_hrs[:, 1], \
                linewidth=1, linestyle='--', marker='>', markersize=markersize, \
                markevery=markerskip, color='k', label='HRS')
        #-----------------------------------------
        ax.plot(seamount_vrs[:, 0], seamount_vrs[:, 1], \
                linewidth=1, linestyle='--', marker='s', markersize=markersize, \
                markevery=markerskip, color='k', label='VRS')
        #-----------------------------------------
        ax.plot(seamount_kmke[:, 0], seamount_kmke[:, 1], \
                linewidth=2, linestyle='-', marker='o', markersize=markersize, \
                markevery=markerskip, color='r', label='KmKe')  
        #-----------------------------------------
        ax.plot(seamount_vebf[:, 0], seamount_vebf[:, 1], \
                linewidth=2, linestyle='-', marker='x', markersize=markersize, \
                markevery=markerskip, color='b', label='VEBF')
        #-----------------------------------------        
        ax.set_xlabel("z",fontsize=18)
        ax.set_ylabel("KmKe",fontsize=18)
        ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)

        ax.set_xlim([0, -5000])
        ax.set_ylim([min_KmKe_ylim, max_KmKe_ylim])

        xticks= np.array([-5000, -4000, -3000, -2000, -1000,])
        # xticklabels = np.array([r'$5000$', r'$4000$', r'$3000$', r'$2000$', r'$1000$'])
        xticklabels = np.array(['5000', '4000', '3000', '2000', '1000'])
        ax.set_xticks(xticks) 
        ax.set_xticklabels(xticklabels)

        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0)) 
        plt.tight_layout()
        xleft, xright = ax.get_xlim()
        ybottom, ytop = ax.get_ylim()
        ratio = 0.5
        ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
        ax.legend(loc='best',  prop={'size': 10})
        fig.savefig(str(indxRange) +'/seamount_avg_kmke_decomp_' + suffix + '.pdf')
        # fig.savefig(str(indxRange) + '/seamount_avg_kmke_decomp_' + suffix + '.png')
        plt.close() # clear figures to avoid overwriting in loop
        #############################     

