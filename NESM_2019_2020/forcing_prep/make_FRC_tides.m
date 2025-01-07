% build frc file for GOM_hycom run (only tides will be included)
clear;
Yorig = 2018;
L=299; M=277;
Lp=L+1;
Mp=M+1;

dir_path = '/storage/home/hcoda1/9/paghor3/scratch/work/NESM_2019_2020/'
 
grdname=[dir_path, 'forcing/','NESM_grd.nc'];

frcout = [dir_path, 'forcing/NESM_frc_tides_20190101_20200229.nc'];
ncid=netcdf.create(frcout,'NETCDF4');
xi_u=netcdf.defDim(ncid,'xi_u',L);
xi_v=netcdf.defDim(ncid,'xi_v',Lp);
eta_u=netcdf.defDim(ncid,'eta_u',Mp);
eta_v=netcdf.defDim(ncid,'eta_v',M);
xi_rho=netcdf.defDim(ncid,'xi_rho',Lp);
eta_rho=netcdf.defDim(ncid,'eta_rho',Mp);
netcdf.endDef(ncid);
netcdf.close(ncid)

% read model grid
lonr=ncread(grdname,'lon_rho'); latr=ncread(grdname,'lat_rho');
lonu=ncread(grdname,'lon_u'); latu=ncread(grdname,'lat_u');
lonv=ncread(grdname,'lon_v'); latv=ncread(grdname,'lat_v');
maskr=ncread(grdname,'mask_rho'); masku=ncread(grdname,'mask_u');
maskv=ncread(grdname,'mask_v');

% add tide forcing
disp('Adding tides ...')
% addpath('~/croco/croco_tools')
% addpath([dir_path, 'croco_tools'])
addpath('/storage/home/hcoda1/9/paghor3/croco/croco_tools/')
start
make_tides
