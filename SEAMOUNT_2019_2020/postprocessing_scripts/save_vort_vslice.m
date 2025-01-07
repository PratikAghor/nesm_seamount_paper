%
% save horizontal slices of vorticity for selected vertical levles (vlevels)
% and time-stamps
%
%-------------------------------------------------------
start_paths
Omega = 7.2921e-5; % rotation rate of Earth s^-1
%-------------------------------------------------------
%%Loop to calculate averages if need (modify as needed)
indxRange = 469:1691; % What time indices do you need?
init = 0; %Index of first time period
% Aghor's script
[~, Nt] = size(indxRange);
hke_arr = zeros(Nt, NumLayers); % get avg w as a function of height and time
time_arr = zeros(Nt, 1);

%% plot options
legendCell = [];
lat_arr = ([lat_rho_vec(90, 1), lat_rho_vec(112, 1), lat_rho_vec(135, 1)]);
lonsec = lon_rho_vec;

lon_arr = ([lon_rho_vec(1, 120), lon_rho_vec(1, 140), lon_rho_vec(1, 165), lon_rho_vec(1, 180)]); % const lon values to 
latsec = lat_rho_vec;

Nz = NumLayers; % no. of vertical levels
vort_t_vslice = zeros(Nt, length(lat_arr), Nx, Nz);

vort_t_vslice_const_lon = zeros(Nt, Ny, length(lon_arr), Nz);

%%
for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    time = ncread(fname, 'time');
    time_arr(nt) = time;
    vort_t_slices = zeros(Ny, Nx, NumLayers);

     % Calculate vorticity for each layer
    for vlevel = 1:NumLayers
	[lat_rho, lon, mask, vort_t_slices(:, :, vlevel)]=get_vort(hisfile,gridfile,tindex,vlevel,coef);
    end

    for k = 1:length(lat_arr) % index for layer numbers
        f0 = 2*Omega*sin(lat_arr(k)); % Coriolis freq.
        lat_val = lat_arr(k); % actual const latitude
        lat_idx = find(abs(lat_rho_vec-lat_arr(k))<1e-3); % find idx of lat in lat_rho_vec
        % Normalize with coriolis freq.
        vort_t_vslice(nt, k, :, :) = squeeze(vort_t_slices(lat_idx, :, :))./f0;   
    end

    for j =1:length(lon_arr)
	f0 = 2*Omega*sin(lat_rho_vec(Ny/2, 1)); % take average f0 over latitudes
	lon_val = lon_arr(j); % actual const longitude
	lon_idx = 172; 
	vort_t_vslice_const_lon(nt, :, j, :) = squeeze(vort_t_slices(:, lon_idx, :))./f0;
    end
end

%%
% const lat vslice
% save data into a netcdf file
vort_t_filename = strcat('seamount_2019_2020_vort_t_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
ncid = netcdf.create(vort_t_filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
y_len = netcdf.defDim(ncid, 'Nlat', length(lat_arr));
x_len = netcdf.defDim(ncid, 'Nx', Nx);
depth_len = netcdf.defDim(ncid, 'Nz', Nz);

netcdf.close(ncid);
% define variables and attributes
ncid = netcdf.open(vort_t_filename,'WRITE');
netcdf.reDef(ncid);

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');

lat_varid = netcdf.defVar(ncid, 'lat', 'double', [y_len]);
netcdf.putAtt(ncid, lat_varid, 'description', 'const lat values where vertical slice is taken');
netcdf.putAtt(ncid, lat_varid, 'units', '--');
netcdf.putAtt(ncid, lat_varid, 'array dimensions', size(y_len));

vort_varid = netcdf.defVar(ncid, 'vort_t_vslice', 'double', [t_len y_len x_len depth_len]);
netcdf.putAtt(ncid, vort_varid, 'description', 'time series of vertical slices of vorticity at fixed lat values');
netcdf.putAtt(ncid, vort_varid, 'units', 'nondimensionalized by Corilis freq.');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(vort_t_vslice));

% close define mode
netcdf.endDef(ncid);

% put values
netcdf.putVar(ncid, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
netcdf.putVar(ncid, lat_varid, lat_arr);
netcdf.putVar(ncid, vort_varid, vort_t_vslice);
% close netcdf file
netcdf.close(ncid);
%%
%% const lon vort_vslice
% save data into a netcdf file
vort_t_filename = strcat('seamount_2019_2020_vort_t_vslice_const_lon_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
ncid2 = netcdf.create(vort_t_filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid2, 'Nt', Nt);
y_len = netcdf.defDim(ncid2, 'Ny', Ny);
x_len = netcdf.defDim(ncid2, 'Nlon', length(lon_arr));
depth_len = netcdf.defDim(ncid2, 'Nz', Nz);

netcdf.close(ncid2);
% define variables and attributes
ncid2 = netcdf.open(vort_t_filename,'WRITE');
netcdf.reDef(ncid2);

time_varid = netcdf.defVar(ncid2, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid2, time_varid, 'description', 'time array');
netcdf.putAtt(ncid2, time_varid, 'units', 's');

lon_varid = netcdf.defVar(ncid2, 'lon', 'double', [x_len]);
netcdf.putAtt(ncid2, lat_varid, 'description', 'const lon values where vertical slice is taken');
netcdf.putAtt(ncid2, lat_varid, 'units', '--');
netcdf.putAtt(ncid, lat_varid, 'array dimensions', size(x_len));

vort_varid = netcdf.defVar(ncid2, 'vort_t_vslice', 'double', [t_len y_len x_len depth_len]);
netcdf.putAtt(ncid2, vort_varid, 'description', 'time series of vertical slices of vorticity at fixed lon values');
netcdf.putAtt(ncid2, vort_varid, 'units', 'nondimensionalized by Corilis freq. (at latitude corr. to mid-box)');
netcdf.putAtt(ncid2, vort_varid, 'array dimensions', size(vort_t_vslice_const_lon));

% close define mode
netcdf.endDef(ncid2);

% put values
netcdf.putVar(ncid2, time_varid, time_arr); % start from 0
netcdf.putVar(ncid2, lon_varid, lon_arr);
netcdf.putVar(ncid2, vort_varid, vort_t_vslice_const_lon);
% close netcdf file
netcdf.close(ncid2);
%%


disp('Done creating vort_t for selected depth layers!')

%%

   
