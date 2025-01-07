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
at z = -1000 and plot spectra
'''

#############################
suffix = 'velo'
#############################
vlevel = -4500
filename1="seamount_2019_2020_speed_t_hslice_vlevel_" + str(vlevel) + "_nt_0_3382.nc"
# filename2="seamount_2019_speed_t_hslice_vlevel_0_nt_1000_1999.nc"
# filename3="seamount_2019_speed_t_hslice_vlevel_0_nt_2000_2902.nc"

data1 = Dataset(filename1, mode="r")
# data2 = Dataset(filename2, mode="r")
# data3 = Dataset(filename3, mode="r")

print("len(data.variables.keys()) = ", len(data1.variables.keys()))
# print("(data.variables.keys()) = ", (data1.variables.keys()))

time_data1 = data1.variables['t_arr']
# vlevel = data1.variables['vlevel']

# this was a typo while saving data in , used 'u_t' instead of 'speed_t' in save_speed_hslice.m
speed_t_1 = data1.variables['u_t'] 

# print("vlevels = ", vlevel[:])
print("vlevels = ", vlevel)

# time_data2 = data2.variables['t_arr']
# speed_t_2 = data2.variables['u_t']

# time_data3 = data3.variables['t_arr']
# speed_t_3 = data3.variables['u_t']

[Nz, Nx, Ny, Nt_file] = np.shape(speed_t_1)
print("np.shape(velo_t_1) = (Nz, Nx, Ny, Nt_file) = ", Nz, ", ", Nx, ", ",  Ny, ", ", Nt_file)


# lon_idx = 158
# lat_idx = 112 
dtsave = time_data1[1] - time_data1[0]
print("dtsave", dtsave)

#############################
#############################
time_data = time_data1
# time_data = np.concatenate((time_data, time_data2[:]), axis=0)
# time_data = np.concatenate((time_data, time_data3[:]), axis=0)

time_data = time_data - time_data[0]

# print("len(timeseries_data) = ", len(timeseries_data))
print("len(time_data) = ", len(time_data))

Nt = len(time_data)
nsave = dtsave

n_files_per_day = 8; # data was saved every 3 hours, 8 savefiles per day
Nt_avg = int(Nt/n_files_per_day)

print("Nt = ", Nt, "n_files_per_day = ", n_files_per_day, \
        "Nt_avg = int(Nt/n_files_per_day) =", Nt_avg)


velo_xymean = zeros((Nz, Nt))
velo_xymean_full = zeros((Nz, Nt))

velo_avg_full = zeros((Nz, Nt_avg)) # declare avg_w, xyt mean
velo_avg = zeros((Nz, Nt_avg)) # declare avg_w, xyt mean
t_avg = zeros((Nt_avg)) # declare avg_w, xyt mean

lon_rho_idx_min = 120
lon_rho_idx_max = 195
lat_rho_idx_min = 55
lat_rho_idx_max = 155

# weight_xyavg = ((lat_rho_idx_max - lat_rho_idx_min)*(lon_rho_idx_max - lon_rho_idx_min))/(Ny*Nx)

velo_hslice = zeros((Ny, Nx, Nt))
velo_hslice_box = zeros((lat_rho_idx_max - lat_rho_idx_min,lon_rho_idx_max - lon_rho_idx_min, Nt))

#############################
# idx i goes in y dir -> lat
# idx j goes in x dir -> lon

for k in range(0, Nz):
        # only pick up w from around Atlantis II
        velo_hslice = speed_t_1[k, :, :, :]
        # velo_hslice = np.concatenate((velo_hslice, speed_t_2[k, :, :, :]), axis=2)
        # velo_hslice = np.concatenate((velo_hslice, speed_t_3[k, :, :, :]), axis=2)
        
        velo_hslice_box = velo_hslice[lon_rho_idx_min:lon_rho_idx_max, \
        lat_rho_idx_min:lat_rho_idx_max, :]
        # velo_hslice_box = velo_hslice
        
        print("np.shape(velo_hslice) = ",  np.shape(velo_hslice))
        print("np.shape(velo_hslice_box) = ",  np.shape(velo_hslice_box))

        # xy avg
        for l in range(0, Nt):
                velo_xymean_full[k, l] = (nanmean(velo_hslice[:, :, l])) 
                velo_xymean[k, l] = (nanmean(velo_hslice_box[:, :, l])) 

        # time avg
        for ll in range(0, Nt_avg):
                velo_avg_full[k, ll] =  (nanmean(velo_xymean_full[k, ll*n_files_per_day:(ll+1)*(n_files_per_day)]))

                velo_avg[k, ll] =  (nanmean(velo_xymean[k, ll*n_files_per_day:(ll+1)*(n_files_per_day)]))
                t_avg[ll] =  time_data[ll*n_files_per_day]
                # print("w_avg[k, ll] = ", w_avg[k, ll])      

np.savetxt("../../compare/avg_surf_velo/seamount_avg_full_surf_velo_vlevel_" \
        + str(vlevel) + ".asc", velo_avg_full[~np.isnan(velo_avg_full)])
np.savetxt("../../compare/avg_surf_velo/seamount_avg_surf_velo_vlevel_" + str(vlevel) + ".asc", velo_avg[~np.isnan(velo_avg)])
print("done creating xyt averages, plotting now")
#############################
