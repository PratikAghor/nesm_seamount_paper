% read and plot instantaneous vorticity hslice (2d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

clear all
close all
start_paths
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
nt = 2045;
vlevel = -3000;

filename = strcat('seamount_2019_2020_vort_2d_hslice_nt_',string(nt), '_vlevel_', string(vlevel), '.nc');

vort_2d = ncread(filename, 'vort_2d');
%--------------------------------------------------------------------------
% only to get the dates in the title!
% vort_vlevel = -4500;
% vort_indxRange = 0:3382; 
% [~, vort_Nt] = size(vort_indxRange);
% nt0=vort_indxRange(1);
% 
% vort_t_filename = strcat('../seamount_2019_2020_vort_t_hslice_vlevel_', string(vort_vlevel(1)), '_nt_',string(vort_indxRange(1)), '_', ...
%     string(vort_indxRange(vort_Nt)), '.nc');
% t_arr = ncread(vort_t_filename, 't_arr');
% t_i = t_arr(indxRange(1) - nt0);
% t_e = t_arr(indxRange(Nt) - nt0);
% date_i = char(string(datetime(2017, 12, 31, 23, 29, t_i)));
% date_e = char(string(datetime(2017, 12, 31, 23, 29, t_e)));
% 
% date_i = date_i(1:6);
% date_e = date_e(1:11);
% 
% date = strcat(date_i, '--', date_e);
%-------------------------------------------------------


%--------------------------------------------------------------------------
%% plot vort hslice
figure1 = figure(3);
% [latlim, lonlim] = geoquadline(lat_rho, lon_rho); geoquadline requires
% Mapping toolbox
% Create axes
ax1 = axes('Parent', figure1, 'YMinorTick','on',...
    'LineWidth',3,...
    'FontSize',24);

m_proj('miller','long', lonlim,'lat', latlim);


hold on;
zMin = -1; % min(min(vort)); 
zMax = -zMin;
h1 = m_image(lon_rho_vec, lat_rho_vec, vort_2d);
colormap(cmocean('balance')); colorbar; 
clim([zMin zMax]);
% colormap(ax1,b2r(zMin,zMax));  colorbar;
% colormap(ax1, "jet"); clim([zMin zMax]); colorbar;
% colormap(ax1, whitejet); clim([zMin zMax]); colorbar;
% freezeColors; hold on;

[h2, c2] = m_contourf(lon_rho, lat_rho, -depth, [vlevel vlevel], ...
    'LineWidth', 3, 'EdgeColor', [0 0 0], 'FaceColor', [1  1 1], 'FaceAlpha', 1);
% c2.FaceColor = [1 1 1]; c2.FaceAlpha = 0.3; % opacity
% c2.EdgeColor = [0 0 0];

m_grid('tickdir','in', ...
'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
'ytick',([38 39]), ... % latitude        
'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;

% title(ax1, string(datetime(2017,12,31,0,0,t_arr(vort_t_nt))));
% ax1.TitleHorizontalAlignment = 'center'; % left makes it come to center
hold on;
% title(ax1, string(date));
%%%
set(figure1, 'Visible', 'off'); % stop pop-ups
figname  = [plots_path, 'vort_plots/hslice/seamount_2019_2020_vort_2d_hslice'];

figname = strcat(figname, '_vlevel_', string(vlevel));
figname = strcat(figname, '_nt_', string(nt));
% 
exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above

% exportgraphics(figure1,strcat(figname, '.eps'))
close all;
%--------------------------------------------------------------------------