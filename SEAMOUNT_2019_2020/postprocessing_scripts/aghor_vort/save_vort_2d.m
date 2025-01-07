% save 2d hslice of vorticity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all
close all
start_paths
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
% read files and calculate mean (bar) values of u, v
nt=3045;

vlevel_arr = ([-3000]);
vlevel=vlevel_arr(1);

vort_3d_filename = strcat('seamount_2019_2020_instant_vort_3d_nt_', string(nt), '.nc');
vort_3d = ncread(vort_3d_filename, 'vort_3d');
vort_3d_sigma = vort_3d./f0;
vort_2d = zeros(Ny, Nx);
%----------------------------------------------------------------------------------

%
% Verticaly integrate
%
time=0;
fname = sprintf([dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], time);
hisfile = fname;
zeta = pagetranspose(ncread(hisfile, 'zeta'));

zr=zlevs(depth,zeta,theta_s,theta_b,hc,N,'r',Vtransform);
zw=zlevs(depth,zeta,theta_s,theta_b,hc,N,'w',Vtransform);
% mask=maskr;
% mask(mask==0)=NaN;
% [HRS2D,h0]=vintegr2(HRS,zw,zr,z01,z02);
%----------------------------------------------------------------------------------
% take a horizontal slice at a given vlevel, modified from get_hslice.m
%
% Read the 3d matrix and do the interpolation
%

vort_2d = vinterp(vort_3d_sigma,zr,vlevel);
%----------------------------------------------------------------------------------

%
% Save time average
%
% save data into a netcdf file
filename = strcat('seamount_2019_2020_vort_2d_hslice_nt_',string(nt), '_vlevel_', string(vlevel), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
% t_len = netcdf.defDim(ncid, 'Nt', Nt);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

vebf_varid = netcdf.defVar(ncid, 'vort_2d', 'double', [y_len x_len]);
netcdf.putAtt(ncid, vebf_varid, 'description', 'time avg of hslice of vort at a given vlevel');
netcdf.putAtt(ncid, vebf_varid, 'units', 'dimless');
netcdf.putAtt(ncid, vebf_varid, 'array dimensions', size(vort_2d));

vlevel_varid = netcdf.defVar(ncid, 'vlevel', 'double', one);
netcdf.putAtt(ncid, vlevel_varid, 'description', 'vlevel for horizontal slice (<0 means depth)');
netcdf.putAtt(ncid, vlevel_varid, 'units', 'm');
netcdf.putAtt(ncid, vlevel_varid, 'array dimensions', size(z01));

% z02_varid = netcdf.defVar(ncid, 'z02', 'double', one);
% netcdf.putAtt(ncid, z02_varid, 'description', 'z02 used for vertical integration');
% netcdf.putAtt(ncid, z02_varid, 'units', 'm');
% netcdf.putAtt(ncid, z02_varid, 'array dimensions', size(z02));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, vlevel_varid, vlevel); % start from 0
% netcdf.putVar(ncid, z02_varid, z02); % start from 0
netcdf.putVar(ncid, vebf_varid, vort_2d); % start from 0
% close netcdf file
netcdf.close(ncid);
disp('Done creating vort_2d at a given vlevel!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


