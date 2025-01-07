%
% save 2d vertical slices of instantaneous 3d vorticity field
%
%-------------------------------------------------------
start_paths
Omega = 7.2921e-5; % rotation rate of Earth s^-1
%-------------------------------------------------------
%%Loop to calculate averages if need (modify as needed)
nt=2132;
% Aghor's script
%--------------------------------------------------------
vort_3d_filename = strcat('seamount_2019_2020_instant_vort_3d_nt_', string(nt), '.nc');
vort_3d = ncread(vort_3d_filename, 'vort_3d');
%--------------------------------------------------------
%% plot options
legendCell = [];
lat_idx_arr = ([56, 79, 90, 101]);
lat_arr = zeros(length(lat_idx_arr), 1); % const lat values to save vslice at
for i = 1:length(lat_idx_arr)
	lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
lat_arr
lonsec = lon_rho_vec;
%--------------------------------------------------------
% lon_idx_arr = ([150, 158, 166, 172]);
% lon_arr = zeros(length(lon_idx_arr), 1);
% for j = 1:length(lon_idx_arr)
% 	lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j)); % const lon values to save vslice at
% end
% latsec = lat_rho_vec;
%--------------------------------------------------------
Nz = NumLayers; % no. of vertical levels
vort_t_vslice = zeros(Nz, length(lat_arr), Nx);

% pv_t_vslice_const_lon = zeros(Nt, Ny, length(lon_arr), Nz);
%--------------------------------------------------------
%%
for i = 1:length(lat_arr) % index for layer numbers
    f0 = 2*Omega*sin(lat_arr(i)); % Coriolis freq.
    lat_val = lat_arr(i); % actual const latitude
    % lat_idx = find(abs(lat_rho_vec-lat_arr(i))<1e-3); % find idx of lat in lat_rho_vec
    lat_idx = lat_idx_arr(i)
% Normalize with coriolis freq.
    vort_t_vslice(:, i, :) = squeeze(vort_3d(:, lat_idx, :))./f0;   
end

    % for j =1:length(lon_arr)
	% f0 = 2*Omega*sin(lat_rho_vec(Ny/2, 1)); % take average f0 over latitudes
	% % lon_val = lon_arr(j); % actual const longitude
	% lon_idx = lon_idx_arr(j); 
	% pv_t_vslice_const_lon(nt, :, j, :) = squeeze(pv_3d(:, lon_idx, :));
    % end
%%
% const lat vslice
% save data into a netcdf file
vort_t_filename = strcat('seamount_2019_2020_instant_vort_vslice_const_lat_nt_',string(nt),'.nc');
ncid = netcdf.create(vort_t_filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
y_len = netcdf.defDim(ncid, 'Nlat', length(lat_arr));
x_len = netcdf.defDim(ncid, 'Nx', Nx);
depth_len = netcdf.defDim(ncid, 'Nz', Nz);

netcdf.close(ncid);
% define variables and attributes
ncid = netcdf.open(vort_t_filename,'WRITE');
netcdf.reDef(ncid);

% time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
% netcdf.putAtt(ncid, time_varid, 'description', 'time array');
% netcdf.putAtt(ncid, time_varid, 'units', 's');

lat_idx_varid = netcdf.defVar(ncid, 'lat_idx', 'double', [y_len]);
netcdf.putAtt(ncid, lat_idx_varid, 'description', 'indices in lat_rho_vec for const lat values where vertical epv slice is taken');
netcdf.putAtt(ncid, lat_idx_varid, 'units', '--');
netcdf.putAtt(ncid, lat_idx_varid, 'array dimensions', size(y_len));

lat_varid = netcdf.defVar(ncid, 'lat', 'double', [y_len]);
netcdf.putAtt(ncid, lat_varid, 'description', 'const lat values where vertical epv slice is taken');
netcdf.putAtt(ncid, lat_varid, 'units', '--');
netcdf.putAtt(ncid, lat_varid, 'array dimensions', size(y_len));

vort_varid = netcdf.defVar(ncid, 'vort_t_vslice', 'double', [depth_len y_len x_len]);
netcdf.putAtt(ncid, vort_varid, 'description', 'normalized vertical slice of relative vorticity at fixed lat values');
netcdf.putAtt(ncid, vort_varid, 'units', 'nondimensionalized by Corilis freq.');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(vort_t_vslice));

% close define mode
netcdf.endDef(ncid);

% put values
% netcdf.putVar(ncid, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
netcdf.putVar(ncid, lat_idx_varid, lat_idx_arr);
netcdf.putVar(ncid, lat_varid, lat_arr);
netcdf.putVar(ncid, vort_varid, vort_t_vslice);
% close netcdf file
netcdf.close(ncid);
%%
%% const lon vort_vslice
% save data into a netcdf file
% vort_t_filename = strcat('seamount_2019_2020_pv_t_vslice_const_lon_nt_',string(indxRange(1)), '_', ...
%     string(indxRange(Nt)),'.nc');
% ncid2 = netcdf.create(vort_t_filename,'CLOBBER');
% % freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
% t_len = netcdf.defDim(ncid2, 'Nt', Nt);
% y_len = netcdf.defDim(ncid2, 'Ny', Ny);
% x_len = netcdf.defDim(ncid2, 'Nlon', length(lon_arr));
% depth_len = netcdf.defDim(ncid2, 'Nz', Nz);
% 
% netcdf.close(ncid2);
% % define variables and attributes
% ncid2 = netcdf.open(vort_t_filename,'WRITE');
% netcdf.reDef(ncid2);
% 
% time_varid = netcdf.defVar(ncid2, 't_arr', 'double', [t_len]);
% netcdf.putAtt(ncid2, time_varid, 'description', 'time array');
% netcdf.putAtt(ncid2, time_varid, 'units', 's');
% 
% 
% lon_idx_varid = netcdf.defVar(ncid2, 'lon_idx', 'double', [x_len]);
% netcdf.putAtt(ncid2, lon_idx_varid, 'description', 'indices in lon_rho_vct for const lon values where vertical epv slice is taken');
% netcdf.putAtt(ncid2, lon_idx_varid, 'units', '--');
% netcdf.putAtt(ncid, lon_idx_varid, 'array dimensions', size(x_len));
% 
% lon_varid = netcdf.defVar(ncid2, 'lon', 'double', [x_len]);
% netcdf.putAtt(ncid2, lon_varid, 'description', 'const lon values where vertical epv slice is taken');
% netcdf.putAtt(ncid2, lon_varid, 'units', '--');
% netcdf.putAtt(ncid, lon_varid, 'array dimensions', size(x_len));
% 
% vort_varid = netcdf.defVar(ncid2, 'pv_t_vslice', 'double', [t_len y_len x_len depth_len]);
% netcdf.putAtt(ncid2, vort_varid, 'description', 'time series of normalized vertical slices of epv at fixed lon values');
% netcdf.putAtt(ncid2, vort_varid, 'units', 'nondimensionalized by Corilis freq. (at latitude corr. to mid-box)');
% netcdf.putAtt(ncid2, vort_varid, 'array dimensions', size(pv_t_vslice_const_lon));
% 
% % close define mode
% netcdf.endDef(ncid2);
% 
% % put values
% netcdf.putVar(ncid2, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid2, lon_idx_varid, lon_idx_arr);
% netcdf.putVar(ncid2, lon_varid, lon_arr);
% netcdf.putVar(ncid2, vort_varid, pv_t_vslice_const_lon);
% % close netcdf file
% netcdf.close(ncid2);
% %%


disp('Done creating vort_t vslices!')

%%

   
