%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Edit a CROCO grid file from NESM_grd.nc
%  to replace NESM by a single ideal seamount 
%  written by: Pratik Aghor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NESM_grdfile='NESM_grd.nc';
SEAMOUNT_grdfile='SEAMOUNT_grd.nc';

disp(' ')
disp(['Creating the file : ', SEAMOUNT_grdfile])
disp(' ')
%%
ncid = netcdf.open(NESM_grdfile,'NC_NOWRITE');
ncid2 = netcdf.create(SEAMOUNT_grdfile,'CLOBBER');


info = ncinfo("NESM_grd.nc");
nVar = numel(info.Variables); % number of variables in NESM_grd.nc
nDims = numel(info.Dimensions);

dimname_arr = strings(nDims, 1);
dimlen_arr = zeros(nDims, 1, 'int64');
for i = 0:nDims-1
    [dimname, dimlen] = netcdf.inqDim(ncid, i);
    dimname_arr(i+1) = dimname;
    dimlen_arr(i+1) = dimlen;
end
% % create dimensions for ncid2, order matters since we've order dimlen_arr
 xi_rho = netcdf.defDim(ncid2, 'xi_rho', dimlen_arr(1));
 eta_rho = netcdf.defDim(ncid2, 'eta_rho', dimlen_arr(2));
 
 xi_u = netcdf.defDim(ncid2, 'xi_u', dimlen_arr(3));
 eta_u = netcdf.defDim(ncid2, 'eta_u', dimlen_arr(4));

 xi_v = netcdf.defDim(ncid2, 'xi_v', dimlen_arr(5));
 eta_v = netcdf.defDim(ncid2, 'eta_v', dimlen_arr(6));

 xi_psi = netcdf.defDim(ncid2, 'xi_psi', dimlen_arr(7));
 eta_psi = netcdf.defDim(ncid2, 'eta_psi', dimlen_arr(8));

 one = netcdf.defDim(ncid2, 'one', dimlen_arr(9));
 s_rho = netcdf.defDim(ncid2, 's_rho', dimlen_arr(10));
 s_w = netcdf.defDim(ncid2, 's_w', dimlen_arr(11));

netcdf.close(ncid2);
%%
% define variables and attributes
ncid2 = netcdf.open(SEAMOUNT_grdfile,'WRITE');
netcdf.reDef(ncid2);

% varid runs from 0 to nVar - 1
% varid = 0, 'xl'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 0);
xl_varid = netcdf.defVar(ncid2, 'xl', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 0, n_att); 
    attval = netcdf.getAtt(ncid, 0,attname);
    netcdf.copyAtt(ncid, 0, attname, ncid2, xl_varid);
end
% check
% attname2 = netcdf.inqAttName(ncid2, 0, 0)
% attval2 = netcdf.getAtt(ncid2, 0, attname2)
% varid = 1, 'el'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 1);
el_varid = netcdf.defVar(ncid2, 'el', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 1, n_att); 
    attval = netcdf.getAtt(ncid, 1, attname);
    netcdf.copyAtt(ncid, 1, attname, ncid2, el_varid);
end

% varid = 2, 'depthmin'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 2);
depthmin_varid = netcdf.defVar(ncid2, 'depthmin', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 2, n_att); 
    attval = netcdf.getAtt(ncid, 2, attname);
    netcdf.copyAtt(ncid, 2, attname, ncid2, depthmin_varid);
end

% varid = 3, 'depthmax'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 3);
depthmax_varid = netcdf.defVar(ncid2, 'depthmax', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 3, n_att); 
    attval = netcdf.getAtt(ncid, 3, attname);
    netcdf.copyAtt(ncid, 3, attname, ncid2, depthmax_varid);
end

% varid = 4, 'spherical'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 4);
spherical_varid = netcdf.defVar(ncid2, 'spherical', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 4, n_att); 
    attval = netcdf.getAtt(ncid, 4, attname);
    netcdf.copyAtt(ncid, 4, attname, ncid2, spherical_varid);
end


% varid = 5, 'angle'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 5);
angle_varid = netcdf.defVar(ncid2, 'angle', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 5, n_att); 
    attval = netcdf.getAtt(ncid, 5, attname);
    netcdf.copyAtt(ncid, 5, attname, ncid2, angle_varid);
end

% varid = 6, 'h'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 6);
h_varid = netcdf.defVar(ncid2, 'h', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 6, n_att); 
    attval = netcdf.getAtt(ncid, 6, attname);
    netcdf.copyAtt(ncid, 6, attname, ncid2, h_varid);
end


% varid = 7, 'hraw'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 7);
hraw_varid = netcdf.defVar(ncid2, 'hraw', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 7, n_att); 
    attval = netcdf.getAtt(ncid, 7, attname);
    netcdf.copyAtt(ncid, 7, attname, ncid2, hraw_varid);
end

% varid = 8, 'f'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 8);
f_varid = netcdf.defVar(ncid2, 'f', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 8, n_att); 
    attval = netcdf.getAtt(ncid, 8, attname);
    netcdf.copyAtt(ncid, 8, attname, ncid2, f_varid);
end

% varid = 9, 'pm'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 9);
pm_varid = netcdf.defVar(ncid2, 'pm', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 9, n_att); 
    attval = netcdf.getAtt(ncid, 9, attname);
    netcdf.copyAtt(ncid, 9, attname, ncid2, pm_varid);
end

% varid = 10, 'pn'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 10);
pn_varid = netcdf.defVar(ncid2, 'pn', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 10, n_att); 
    attval = netcdf.getAtt(ncid, 10, attname);
    netcdf.copyAtt(ncid, 10, attname, ncid2, pn_varid);
end


% varid = 11, 'dndx'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 11);
dndx_varid = netcdf.defVar(ncid2, 'dndx', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 11, n_att); 
    attval = netcdf.getAtt(ncid, 11, attname);
    netcdf.copyAtt(ncid, 11, attname, ncid2, dndx_varid);
end

% varid = 12, 'dmde'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 12);
dmde_varid = netcdf.defVar(ncid2, 'dmde', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 12, n_att); 
    attval = netcdf.getAtt(ncid, 12, attname);
    netcdf.copyAtt(ncid, 12, attname, ncid2, dmde_varid);
end

% varid = 13, 'lon_rho'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 13);
lon_rho_varid = netcdf.defVar(ncid2, 'lon_rho', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 13, n_att); 
    attval = netcdf.getAtt(ncid, 13, attname);
    netcdf.copyAtt(ncid, 13, attname, ncid2, lon_rho_varid);
end

% varid = 14, 'lon_u'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 14);
lon_u_varid = netcdf.defVar(ncid2, 'lon_u', xtype, [xi_u eta_u]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 14, n_att); 
    attval = netcdf.getAtt(ncid, 14, attname);
    netcdf.copyAtt(ncid, 14, attname, ncid2, lon_u_varid);
end

% varid = 15, 'lon_v'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 15);
lon_v_varid = netcdf.defVar(ncid2, 'lon_v', xtype, [xi_v eta_v]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 15, n_att); 
    attval = netcdf.getAtt(ncid, 15, attname);
    netcdf.copyAtt(ncid, 15, attname, ncid2, lon_v_varid);
end

% varid = 16, 'lon_psi'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 16);
lon_psi_varid = netcdf.defVar(ncid2, 'lon_psi', xtype, [xi_psi eta_psi]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 16, n_att); 
    attval = netcdf.getAtt(ncid, 16, attname);
    netcdf.copyAtt(ncid, 16, attname, ncid2, lon_psi_varid);
end

% varid = 17, 'lat_rho'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 17);
lat_rho_varid = netcdf.defVar(ncid2, 'lat_rho', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 17, n_att); 
    attval = netcdf.getAtt(ncid, 17, attname);
    netcdf.copyAtt(ncid, 17, attname, ncid2, lat_rho_varid);
end

% varid = 18, 'lat_u'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 18);
lat_u_varid = netcdf.defVar(ncid2, 'lat_u', xtype, [xi_u eta_u]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 18, n_att); 
    attval = netcdf.getAtt(ncid, 18, attname);
    netcdf.copyAtt(ncid, 18, attname, ncid2, lat_u_varid);
end

% varid = 19, 'lat_v'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 19);
lat_v_varid = netcdf.defVar(ncid2, 'lat_v', xtype, [xi_v eta_v]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 19, n_att); 
    attval = netcdf.getAtt(ncid, 19, attname);
    netcdf.copyAtt(ncid, 19, attname, ncid2, lat_v_varid);
end

% varid = 20, 'lat_psi'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 20);
lat_psi_varid = netcdf.defVar(ncid2, 'lat_psi', xtype, [xi_psi eta_psi]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 20, n_att); 
    attval = netcdf.getAtt(ncid, 20, attname);
    netcdf.copyAtt(ncid, 20, attname, ncid2, lat_psi_varid);
end
% varid = 21, 'mask_rho'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 21);
mask_rho_varid = netcdf.defVar(ncid2, 'mask_rho', xtype, [xi_rho eta_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 21, n_att); 
    attval = netcdf.getAtt(ncid, 21, attname);
    netcdf.copyAtt(ncid, 21, attname, ncid2, mask_rho_varid);
end

% varid = 22, 'mask_u'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 22);
mask_u_varid = netcdf.defVar(ncid2, 'mask_u', xtype, [xi_u eta_u]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 22, n_att); 
    attval = netcdf.getAtt(ncid, 22, attname);
    netcdf.copyAtt(ncid, 22, attname, ncid2, mask_u_varid);
end

% varid = 23, 'mask_v'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 23);
mask_v_varid = netcdf.defVar(ncid2, 'mask_v', xtype, [xi_v eta_v]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 23, n_att); 
    attval = netcdf.getAtt(ncid, 23, attname);
    netcdf.copyAtt(ncid, 23, attname, ncid2, mask_v_varid);
end

% varid = 24, 'mask_psi'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 24);
mask_psi_varid = netcdf.defVar(ncid2, 'mask_psi', xtype, [xi_psi eta_psi]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 24, n_att); 
    attval = netcdf.getAtt(ncid, 24, attname);
    netcdf.copyAtt(ncid, 24, attname, ncid2, mask_psi_varid);
end

% varid = 25, 'hc'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 25);
hc_varid = netcdf.defVar(ncid2, 'hc', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 25, n_att); 
    attval = netcdf.getAtt(ncid, 25, attname);
    netcdf.copyAtt(ncid, 25, attname, ncid2, hc_varid);
end

% varid = 26, 'Cs_r'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 26);
Cs_r_varid = netcdf.defVar(ncid2, 'Cs_r', xtype, [s_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 26, n_att); 
    attval = netcdf.getAtt(ncid, 26, attname);
    netcdf.copyAtt(ncid, 26, attname, ncid2, Cs_r_varid);
end

% varid = 27, 'Cs_w'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 27);
Cs_w_varid = netcdf.defVar(ncid2, 'Cs_w', xtype, [s_w]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 27, n_att); 
    attval = netcdf.getAtt(ncid, 27, attname);
    netcdf.copyAtt(ncid, 27, attname, ncid2, Cs_w_varid);
end

% varid = 28, 's_rho'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 28);
s_rho_varid = netcdf.defVar(ncid2, 's_rho', xtype, [s_rho]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 28, n_att); 
    attval = netcdf.getAtt(ncid, 28, attname);
    netcdf.copyAtt(ncid, 28, attname, ncid2, s_rho_varid);
end
% varid = 29, 's_w'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 29);
s_w_varid = netcdf.defVar(ncid2, 's_w', xtype, [s_w]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 29, n_att); 
    attval = netcdf.getAtt(ncid, 29, attname);
    netcdf.copyAtt(ncid, 29, attname, ncid2, s_w_varid);
end

% varid = 30, 'Vtransform'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 30);
Vtransform_varid = netcdf.defVar(ncid2, 'Vtransform', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 30, n_att); 
    attval = netcdf.getAtt(ncid, 30, attname);
    netcdf.copyAtt(ncid, 30, attname, ncid2, Vtransform_varid);
end

% varid = 31, 'Vstretching'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 31);
Vstretching_varid = netcdf.defVar(ncid2, 'Vstretching', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 31, n_att); 
    attval = netcdf.getAtt(ncid, 31, attname);
    netcdf.copyAtt(ncid, 31, attname, ncid2, Vstretching_varid);
end

% varid = 32, 'theta_s'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 32);
theta_s_varid = netcdf.defVar(ncid2, 'theta_s', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 32, n_att); 
    attval = netcdf.getAtt(ncid, 32, attname);
    netcdf.copyAtt(ncid, 32, attname, ncid2, theta_s_varid);
end

% varid = 33, 'theta_b'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 33);
theta_b_varid = netcdf.defVar(ncid2, 'theta_b', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 33, n_att); 
    attval = netcdf.getAtt(ncid, 33, attname);
    netcdf.copyAtt(ncid, 33, attname, ncid2, theta_b_varid);
end

% varid = 34, 'Tcline'
[varname, xtype, dimids, atts] = netcdf.inqVar(ncid, 34);
Tcline_varid = netcdf.defVar(ncid2, 'Tcline', xtype, [one]);
% copy attributes
for n_att = 0:atts-1
    attname = netcdf.inqAttName(ncid, 34, n_att); 
    attval = netcdf.getAtt(ncid, 34, attname);
    netcdf.copyAtt(ncid, 34, attname, ncid2, Tcline_varid);
end
%%
%
% Create global attributes
%
v_glob = netcdf.getConstant('GLOBAL');
title='SEAMOUNT_s100'
type='ROMS grid file'
netcdf.putAtt(ncid2, v_glob, 'title', title);
netcdf.putAtt(ncid2, v_glob, 'date', date);
netcdf.putAtt(ncid2, v_glob, 'type', type);
% close define mode
netcdf.endDef(ncid2);
%%
% geVar from ncid 
xl = ncread(NESM_grdfile, 'xl');
el = ncread(NESM_grdfile, 'el');
depthmin = ncread(NESM_grdfile, 'depthmin');
depthmax = ncread(NESM_grdfile, 'depthmax');
spherical = ncread(NESM_grdfile, 'spherical');
angle = ncread(NESM_grdfile, 'angle');
f = ncread(NESM_grdfile, 'f');
pm = ncread(NESM_grdfile, 'pm');
pn = ncread(NESM_grdfile, 'pn');
dndx = ncread(NESM_grdfile, 'dndx');
dmde = ncread(NESM_grdfile, 'dmde');
lon_rho = ncread(NESM_grdfile, 'lon_rho');
lon_u = ncread(NESM_grdfile, 'lon_u');
lon_v = ncread(NESM_grdfile, 'lon_v');
lon_psi = ncread(NESM_grdfile, 'lon_psi');
lat_rho = ncread(NESM_grdfile, 'lat_rho');
lat_u = ncread(NESM_grdfile, 'lat_u');
lat_v = ncread(NESM_grdfile, 'lat_v');
lat_psi = ncread(NESM_grdfile, 'lat_psi');
mask_rho = ncread(NESM_grdfile, 'mask_rho');
mask_u = ncread(NESM_grdfile, 'mask_u');
mask_v = ncread(NESM_grdfile, 'mask_v');
mask_psi = ncread(NESM_grdfile, 'mask_psi');
hc = ncread(NESM_grdfile, 'hc');
Cs_r = ncread(NESM_grdfile, 'Cs_r');
Cs_w = ncread(NESM_grdfile, 'Cs_w');
s_rho = ncread(NESM_grdfile, 's_rho');
s_w = ncread(NESM_grdfile, 's_w');
Vtransform = ncread(NESM_grdfile, 'Vtransform');
Vstretching = ncread(NESM_grdfile, 'Vstretching');
theta_s = ncread(NESM_grdfile, 'theta_s');
theta_b = ncread(NESM_grdfile, 'theta_b');
Tcline = ncread(NESM_grdfile, 'Tcline');

% redefine h and hraw
%
% Grid dimensions:
%
lonmin = -65;   % Minimum longitude [degree east]
lonmax = -61;   % Maximum longitude [degree east]
latmin = 37;   % Minimum latitude  [degree north]
latmax = 43;   % Maximum latitude  [degree north]
%
% Atlantis II lat-lon lims
atl_lonmin = -63.6; % -63.5;   % Minimum longitude of Atlantis II  [degree east]
atl_lonmax = -62.7; % -62.84;   % Maximum longitude of Atlantis II [degree east]
atl_latmin = 38;   % Minimum latitude of Atlantis II [degree north]
atl_latmax = 38.9;   % Maximum latitude of Atlantis II [degree north]

% gaussian seamount params
% centre == mean, sx, sy == variances
atl_loncentre = -63.23; % (atl_lonmin + atl_lonmax)/2;
atl_latcentre = 38.5; % (atl_latmin + atl_latmax)/2;
sx = 0.03; %(abs(atl_lonmin - atl_loncentre)/100);
sy = 0.03; %(abs(atl_latmax - atl_latcentre)/100);

h_nesm = ncread(NESM_grdfile, 'h');
h=h_nesm(185, 156)*ones(size(lat_rho));
hmax = h_nesm(185, 156); % depthmax;
height = hmax - 1645; % depthmax - depthmin; % height of the gaussian seamount in meters
% Atlantis II rises to 1645 m 

% h = -hmax*ones(size(lat_rho));
[Ny, Nx] = size(lat_rho);
for i = 1:Ny
    for j = 1:Nx
        local_lat = lat_rho(1, j);
        local_lon = lon_rho(i, 1);
        % if((local_lat >= atl_latmin && local_lat <= atl_latmax) && (local_lon >= atl_lonmin && local_lon <= atl_lonmax))
        %    sprintf('(i, j), (lat, lon) = (%d, %d), (%0.2f, %0.2f)', i, j, local_lat, local_lon)
        % h(i, j) = hmax - (height)*exp(-(local_lon - atl_loncentre)^2/(sx) -(local_lat - atl_latcentre)^2/(sy)); % gaussian seamount
        % end
	h(i, j) = hmax - (height)*exp(-(local_lon - atl_loncentre)^2/(sx) -(local_lat - atl_latcentre)^2/(sy)); % gaussian seamount
    end
end
%%

hraw = h; % ncread(NESM_grdfile, 'hraw');

% putVar vals to ncid2
% write data
netcdf.putVar(ncid2, xl_varid, xl);
netcdf.putVar(ncid2, el_varid, el);
netcdf.putVar(ncid2, depthmin_varid, depthmin);
netcdf.putVar(ncid2, depthmax_varid, depthmax);
netcdf.putVar(ncid2, spherical_varid, spherical);
netcdf.putVar(ncid2, angle_varid, angle);
netcdf.putVar(ncid2, h_varid, h);
netcdf.putVar(ncid2, hraw_varid, hraw);
netcdf.putVar(ncid2, f_varid, f);
netcdf.putVar(ncid2, pm_varid, pm);
netcdf.putVar(ncid2, pn_varid, pn);
netcdf.putVar(ncid2, dndx_varid, dndx);
netcdf.putVar(ncid2, dmde_varid, dmde);
netcdf.putVar(ncid2, lon_rho_varid, lon_rho);
netcdf.putVar(ncid2, lon_u_varid, lon_u);
netcdf.putVar(ncid2, lon_v_varid, lon_v);
netcdf.putVar(ncid2, lon_psi_varid, lon_psi);
netcdf.putVar(ncid2, lat_rho_varid, lat_rho);
netcdf.putVar(ncid2, lat_u_varid, lat_u);
netcdf.putVar(ncid2, lat_v_varid, lat_v);
netcdf.putVar(ncid2, lat_psi_varid, lat_psi);
netcdf.putVar(ncid2, mask_rho_varid, mask_rho);
netcdf.putVar(ncid2, mask_u_varid, mask_u);
netcdf.putVar(ncid2, mask_v_varid, mask_v);
netcdf.putVar(ncid2, mask_psi_varid, mask_psi);

netcdf.putVar(ncid2, hc_varid, hc);
netcdf.putVar(ncid2, Cs_r_varid, Cs_r);
netcdf.putVar(ncid2, Cs_w_varid, Cs_w);
netcdf.putVar(ncid2, s_rho_varid, s_rho);
netcdf.putVar(ncid2, s_w_varid, s_w);
netcdf.putVar(ncid2, Vtransform_varid, Vtransform);
netcdf.putVar(ncid2, Vstretching_varid, Vstretching);
netcdf.putVar(ncid2, theta_s_varid, theta_s);
netcdf.putVar(ncid2, theta_b_varid, theta_b);
netcdf.putVar(ncid2, Tcline_varid, Tcline);
%%
%%
netcdf.close(ncid);
netcdf.close(ncid2);

disp('done creating SEAMOUNT_grd.nc')
