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
# compare KmKe, HRS, VRS, VEBF at const lat-lon as a function of z
# Author: Pratik Aghor

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
lat_arr = np.array([38.7095])
lon_arr = np.array([-63.2558, -62.9992])

if(indxRange=='469_3382'):
        max_KmKe_ylim = 5e-7
        min_KmKe_ylim = -5e-7
else:
        max_KmKe_ylim = 5e-6
        min_KmKe_ylim = -5e-6


markerskip=2
markersize=5
linewidth=3
#############################
for i in range (0,len(lat_arr)):
        for j in range(0, len(lon_arr)):
                lat_lon = '_lat_' + str(lat_arr[i]) + '_lon_' + str(lon_arr[j])
                suffix = lat_lon + '_nt_' + str(indxRange)

                nesm_hrs = np.loadtxt('hrs/nesm_2019_2020_hrs_z'+ lat_lon + '_nt_'+ str(indxRange) + '.asc')
                seamount_hrs = np.loadtxt('hrs/seamount_2019_2020_hrs_z' + lat_lon + '_nt_'+ str(indxRange) + '.asc')

                nesm_vrs = np.loadtxt('vrs/nesm_2019_2020_vrs_z'+ lat_lon +'_nt_'+ str(indxRange) + '.asc')
                seamount_vrs = np.loadtxt('vrs/seamount_2019_2020_vrs_z' + lat_lon + '_nt_'+ str(indxRange) + '.asc')

                nesm_kmke = np.loadtxt('kmke/nesm_2019_2020_kmke_z'+ lat_lon +'_nt_'+ str(indxRange) + '.asc')
                seamount_kmke = np.loadtxt('kmke/seamount_2019_2020_kmke_z' + lat_lon + '_nt_'+ str(indxRange) + '.asc')

                nesm_vebf = np.loadtxt('vebf/nesm_2019_2020_vebf_z' + lat_lon + '_nt_'+ str(indxRange) + '.asc')
                seamount_vebf = np.loadtxt('vebf/seamount_2019_2020_vebf_z' + lat_lon + '_nt_'+ str(indxRange) + '.asc')

                #############################

                # find first occurnace of NaN to set ylim
                seamount_first_nan_idx = np.where(np.isnan(seamount_vebf[:, 1]))[0][0]
                print("seamount_first_nan_idx", seamount_first_nan_idx)
                # set that value to 0 for suit-putting
                seamount_hrs[seamount_first_nan_idx, 1] = 0
                seamount_vrs[seamount_first_nan_idx, 1] = 0
                seamount_kmke[seamount_first_nan_idx, 1] = 0
                seamount_vebf[seamount_first_nan_idx, 1] = 0
                seamount_ybottom = seamount_vebf[seamount_first_nan_idx, 0]

                nesm_first_nan_idx = np.where(np.isnan(nesm_vebf[:, 1]))[0][0]
                print("nesm_first_nan_idx", nesm_first_nan_idx)
                # set that value to 0 for suit-putting
                nesm_hrs[nesm_first_nan_idx, 1] = 0
                nesm_vrs[nesm_first_nan_idx, 1] = 0
                nesm_kmke[nesm_first_nan_idx, 1] = 0
                nesm_vebf[nesm_first_nan_idx, 1] = 0
                nesm_ybottom = nesm_vebf[nesm_first_nan_idx, 0]

                ytop = -3000
                #############################
                fig = plt.figure('seamount_kmke_decomp_' + str(i) + '_' + str(j))
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
                #        linewidth=linewidth, linestyle='-', marker='o', markersize=markersize, \
                #        markevery=markerskip, color='r', label=r'\textrm{KmKe}')        
                # ax.set_xlabel(r"$z$",fontsize=18)
                # ax.set_ylabel(r"$KmKe$",fontsize=18)


                # without latex labels
                #-----------------------------------------
                ax.plot(seamount_hrs[:, 1], seamount_hrs[:, 0], \
                        linewidth=2, linestyle='--', marker='>', markersize=markersize, \
                        markevery=markerskip, color='k', label='HRS')
                #-----------------------------------------
                ax.plot(seamount_vrs[:, 1], seamount_vrs[:, 0], \
                        linewidth=2, linestyle='--', marker='s', markersize=markersize, \
                        markevery=markerskip, color='k', label='VRS')
                #-----------------------------------------
                ax.plot(seamount_kmke[:, 1], seamount_kmke[:, 0], \
                        linewidth=linewidth, linestyle='-', marker='o', markersize=markersize, \
                        markevery=markerskip, color='r', label='KmKe')  
                #-----------------------------------------
                ax.plot(seamount_vebf[:, 1], seamount_vebf[:, 0], \
                        linewidth=linewidth, linestyle='-', marker='x', markersize=markersize, \
                        markevery=markerskip, color='b', label='VEBF')
                #-----------------------------------------        
                ax.set_ylabel("z",fontsize=18)
                # ax.set_xlabel("KmKe, VEBF",fontsize=18)
                ax.axvline(x=0, color='k', linestyle='--', linewidth=0.5)

                ax.set_ylim([0, -5000])
                ax.set_xlim([min_KmKe_ylim, max_KmKe_ylim])

                yticks= np.array([-5000, -4000, -3000, -2000, -1000,])
                yticklabels = np.array([r'$5000$', r'$4000$', r'$3000$', r'$2000$', r'$1000$'])
                # yticklabels = np.array(['5000', '4000', '3000', '2000', '1000'])
                xticks= np.array([min_KmKe_ylim, 0, max_KmKe_ylim])
                ax.set_yticks(yticks) 
                ax.set_yticklabels(yticklabels)
                ax.set_xticks(xticks)

                # plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0)) 
                plt.tight_layout()
                xleft, xright = ax.get_xlim()
                # ybottom, ytop = ax.get_ylim()
                # ax.set_ylim(ax.get_ylim()[::-1])
                
                ax.set_ylim([seamount_ybottom, ytop])
                ratio = 2
                ax.set_aspect(abs((xright-xleft)/(seamount_ybottom-ytop))*ratio)
                ax.legend(loc=2,  prop={'size': 10})
                fig.savefig(str(indxRange) +'/seamount_avg_kmke_decomp' + suffix + '.pdf')
                # fig.savefig(str(indxRange) + '/seamount_avg_kmke_decomp_' + suffix + '.png')
                plt.close() # clear figures to avoid overwriting in loop
                #############################
                #############################
                fig = plt.figure('nesm_kmke_decomp_' + str(i) + '_' + str(j))
                ax = fig.gca()  # Get current axes
                #-----------------------------------------
                # ax.plot(nesm_hrs[:, 0], nesm_hrs[:, 1], \
                #        linewidth=1, linestyle='--', marker='>', markersize=markersize, \
                #        markevery=markerskip, color='k', label=r'\textrm{HRS}')
                #-----------------------------------------
                #ax.plot(nesm_vrs[:, 0], nesm_vrs[:, 1], \
                #        linewidth=1, linestyle='--', marker='s', markersize=markersize, \
                #        markevery=markerskip, color='k', label=r'\textrm{VRS}')
                #-----------------------------------------
                #ax.plot(nesm_kmke[:, 0], nesm_kmke[:, 1], \
                #        linewidth=linewidth, linestyle='-', marker='o', markersize=markersize, \
                #        markevery=markerskip, color='r', label=r'\textrm{KmKe}')        
                # ax.set_xlabel(r"$z$",fontsize=18)
                # ax.set_ylabel(r"$KmKe$",fontsize=18)


                # without latex labels
                #-----------------------------------------
                ax.plot(nesm_hrs[:, 1], nesm_hrs[:, 0], \
                        linewidth=2, linestyle='--', marker='>', markersize=markersize, \
                        markevery=markerskip, color='k', label='HRS')
                #-----------------------------------------
                ax.plot(nesm_vrs[:, 1], nesm_vrs[:, 0], \
                        linewidth=2, linestyle='--', marker='s', markersize=markersize, \
                        markevery=markerskip, color='k', label='VRS')
                #-----------------------------------------
                ax.plot(nesm_kmke[:, 1], nesm_kmke[:, 0], \
                        linewidth=linewidth, linestyle='-', marker='o', markersize=markersize, \
                        markevery=markerskip, color='r', label='KmKe')  
                #-----------------------------------------
                ax.plot(nesm_vebf[:, 1], nesm_vebf[:, 0], \
                        linewidth=linewidth, linestyle='-', marker='x', markersize=markersize, \
                        markevery=markerskip, color='b', label='VEBF')
                #-----------------------------------------        
                ax.set_ylabel("z",fontsize=18)
                # ax.set_xlabel("KmKe, VEBF",fontsize=18)
                # ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
                ax.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
                ax.set_ylim([0, -5000])
                ax.set_xlim([min_KmKe_ylim, max_KmKe_ylim])

                yticks= np.array([-5000, -4000, -3000, -2000, -1000])
                yticklabels = np.array([r'$5000$', r'$4000$', r'$3000$', r'$2000$', r'$1000$'])
                xticks= np.array([min_KmKe_ylim, 0, max_KmKe_ylim])
                # yticklabels = np.array(['5000', '4000', '3000', '2000', '1000'])
                ax.set_yticks(yticks) 
                ax.set_yticklabels(yticklabels)
                ax.set_xticks(xticks) 

                # plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0)) 
                plt.tight_layout()
                xleft, xright = ax.get_xlim()
                # ybottom, ytop = ax.get_ylim()
                # ax.set_ylim(ax.get_ylim()[::-1])
                
                ax.set_ylim([nesm_ybottom, ytop])
                ratio = 2
                ax.set_aspect(abs((xright-xleft)/(nesm_ybottom-ytop))*ratio)
                ax.legend(loc=2,  prop={'size': 10})
                fig.savefig(str(indxRange) +'/nesm_avg_kmke_decomp' + suffix + '.pdf')
                # fig.savefig(str(indxRange) + '/seamount_avg_kmke_decomp_' + suffix + '.png')
                plt.close() # clear figures to avoid overwriting in loop
                #############################