%
% save horizontal slices of w for selected vertical levles (vlevels)
% and time-stamps
%
%-------------------------------------------------------
start_paths
type='r'; % see get_hslice 
%-------------------------------------------------------
%%Loop to calculate averages if need (modify as needed)
indxRange = 0:1690; % What time indices do you need?
init = 0; %Index of first time period
% Aghor's script
[~, Nt] = size(indxRange);
hke_arr = zeros(Nt, NumLayers); % get avg w as a function of height and time
time_arr = zeros(Nt, 1);

%% plot options
% vlevel_skip = 25;
legendCell = [];
% vlevel_arr = ([100, 52, 39, 15, 1]);% NumLayers: -vlevel_skip: 1

% due to size constraints, save only one vlevel at a time
vlevel_arr = ([0]); % vlevel < 0 => depth;
Nz = length(vlevel_arr);
speed_t = zeros(Nt, Ny, Nx, length(vlevel_arr));

% w_t = zeros(Nt, Ny, Nx, length(vlevel_arr));
%%
for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    time = ncread(fname, 'time');
    time_arr(nt) = time;
    for k = 1:Nz % index for layer numbers (since we don't save all layers)
        vlevel = vlevel_arr(k); % actual layer number
	% get uvw at 'r' (rho) points
	[~, ~, ~, speed_t(nt, :, :, k)] = get_speed(hisfile, gridfile, tindex, vlevel, coef);  
    end
    
end

%%
% save data into a netcdf file
% filename = strcat('seamount_2019_speed_t_hslice', '_nt_',string(indxRange(1)), '_', ...
%    string(indxRange(Nt)),'.nc');ncid = netcdf.create(filename,'CLOBBER');
filename = strcat('seamount_2019_2020_speed_t_hslice_vlevel_', string(vlevel_arr(1)) , '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');ncid = netcdf.create(filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
depth_len = netcdf.defDim(ncid, 'Nz', Nz);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');

vlevel_varid = netcdf.defVar(ncid, 'vlevel', 'double', [depth_len]);
netcdf.putAtt(ncid, vlevel_varid, 'description', 'array of vertical layers from the original data');
netcdf.putAtt(ncid, vlevel_varid, 'units', '--');
netcdf.putAtt(ncid, vlevel_varid, 'array dimensions', size(depth_len));

u_varid = netcdf.defVar(ncid, 'u_t', 'double', [t_len y_len x_len depth_len]);
netcdf.putAtt(ncid, u_varid, 'description', 'time series of speed hslice');
netcdf.putAtt(ncid, u_varid, 'units', 'm/s');
netcdf.putAtt(ncid, u_varid, 'array dimensions', size(speed_t));

% w_varid = netcdf.defVar(ncid, 'w_t', 'double', [t_len y_len x_len depth_len]);
% netcdf.putAtt(ncid, w_varid, 'description', 'time series of w hslice');
% netcdf.putAtt(ncid, w_varid, 'units', 'm/s');
% netcdf.putAtt(ncid, w_varid, 'array dimensions', size(w_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
netcdf.putVar(ncid, vlevel_varid, vlevel_arr);
netcdf.putVar(ncid, u_varid, speed_t);
% netcdf.putVar(ncid, w_varid, w_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating speed_t_hslice for selected depth layers!')

%%

   
