%--------------------------------------
% given nt, save 3d vorticity matrix in sigma coordinates 
% Pratik Aghor
%--------------------------------------
start_paths
%--------------------------------------
nt = 3045;

Nz=NumLayers; % no. of vertical levels
vort_3d = zeros(Nz, Ny, Nx);
%----------------------------------------
%%
sprintf(strcat('loading', [dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], 'file'), nt)
fname = sprintf([dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], nt);
hisfile = fname;
% Calculate vorticity for each layer
for vlevel = 1:NumLayers
    [lat_rho, lon, mask, vort_3d(vlevel, :, :)]=get_vort(hisfile,gridfile,tindex,vlevel,coef);
end
%----------------------------------------
% save data into a netcdf file

vort_t_filename = strcat('seamount_2019_2020_instant_vort_3d', '_nt_', string(nt), '.nc');
ncid = netcdf.create(vort_t_filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
depth_len = netcdf.defDim(ncid, 'Nz', NumLayers);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(vort_t_filename,'WRITE');
netcdf.reDef(ncid);

vort_varid = netcdf.defVar(ncid, 'vort_3d', 'double', [depth_len y_len x_len]);
netcdf.putAtt(ncid, vort_varid, 'description', '3d relative vorticity in sigma coords');
netcdf.putAtt(ncid, vort_varid, 'units', 's^-1');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(vort_3d));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, vort_varid, vort_3d);
% close netcdf file
netcdf.close(ncid);
disp('Done creating vort_3d!')
%----------------------------------------------------------------
