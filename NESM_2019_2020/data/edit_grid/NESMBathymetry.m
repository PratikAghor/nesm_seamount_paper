%------------------------------------
% on prometheus
main_dir_path='/home/aghor/';
aghor_extras_path='/home/aghor/aghor/GT/GT_project/aghor_extras/'; % change according to the system
addpath([aghor_extras_path, 'm_map/']);
plots_path = './';
%------------------------------------

% lon = ncread('SEAMOUNT_grd.nc', 'lon_rho');
% lat = ncread('SEAMOUNT_grd.nc', 'lat_rho');
% depth = -ncread('SEAMOUNT_grd.nc', 'h');

lon = pagetranspose(ncread('NESM_grd.nc', 'lon_rho'));
lat = pagetranspose(ncread('NESM_grd.nc', 'lat_rho'));
depth = pagetranspose(ncread('NESM_grd.nc', 'h'));
lon_rho_vec = squeeze(lon(1, :));
lat_rho_vec = squeeze(lat(:, 1));
% 
lon_rho_idx_min = 120;
lon_rho_idx_max = 195;
lat_rho_idx_min = 55;
lat_rho_idx_max = 155;

% %% Plotting
figure1 = figure();
[latlim,lonlim] = geoquadline(lat,lon);
m_proj('miller', 'long', lonlim,'lat', latlim);
zMin = 0;
zMax = 5000;
m_contourf(lon, lat, depth);
clim([zMin, zMax]);
m_grid('tickdir','in', ...
       'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
       'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
       'ytick',([38 39]), ... % latitude        
       'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;

hold on;
m_plot([lon_rho_vec(lon_rho_idx_min) lon_rho_vec(lon_rho_idx_max)],[lat_rho_vec(lat_rho_idx_min) lat_rho_vec(lat_rho_idx_min)], ...
    'r', 'LineWidth', 3);
m_plot([lon_rho_vec(lon_rho_idx_min) lon_rho_vec(lon_rho_idx_max)],[lat_rho_vec(lat_rho_idx_max) lat_rho_vec(lat_rho_idx_max)], ...
    'r', 'LineWidth', 3);

m_plot([lon_rho_vec(lon_rho_idx_min) lon_rho_vec(lon_rho_idx_min)],[lat_rho_vec(lat_rho_idx_min) lat_rho_vec(lat_rho_idx_max)], ...
    'r', 'LineWidth', 3);
m_plot([lon_rho_vec(lon_rho_idx_max) lon_rho_vec(lon_rho_idx_max)],[lat_rho_vec(lat_rho_idx_min) lat_rho_vec(lat_rho_idx_max)], ...
    'r', 'LineWidth', 3);
% aghoredit

% title('NESM Bathymetry')
% title('Atlantis-II bathy')
h = colorbar;
% h.Label.String = "Depth (m)";
% h.Label.Rotation = 270;
% h.Label.VerticalAlignment = "bottom";
h.YTick = [1000 2000 3000 4000 5000];
h.YTickLabel = {'1000', '2000', '3000', '4000', '5000'};
set(figure1, 'Visible', 'off'); % stop pop-ups
figname  = [plots_path, 'nesm_bathy'];
figname = strcat(figname, '.pdf');
exportgraphics(figure1, figname, 'ContentType', 'vector'); % remove extra white space, 2022a and above

%------------------------------------