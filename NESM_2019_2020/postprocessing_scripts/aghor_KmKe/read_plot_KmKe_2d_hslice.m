% 
% Compute the kinetic energy transfer: horizontal Reynolds stress (HRS)
% - Pratik Aghor, modified from get_KmKe.m
% HRS= -[ <up up> dubar/dx + <up vp>  dubar/dy  ....
%          <vp up> dvbar/dx + <vp vp>  dvbar/dy ]
%     =  -<up(up.grad(ubar))+vp(up.grad(vbar))>
%
% Advection operators are used. Much easier in sigma coordinates.
%
% This is the method used in
% Djakouré, S., P. Penven, B. Bourlès, J. Veitch and V. Koné, 
% Coastally trapped eddies in the north of the Gulf of Guinea, 2014, J. Geophys. Res,
% DOI: 10.1002/2014JC010243
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Assume mean and hrs_3d with respect to that mean has been calculated and saved
% only take a horizontal section of the saved hrs_3d at a given depth, 
% and save HRS2D 
 
clear all
close all
start_paths
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
% case_type='instant'
case_type='avg' % 'instant' or 'avg'

vlevel = -4500;

if(strcmp(case_type, 'instant'))
    nt=3040;
    KmKe_2d_file = strcat('nesm_2019_2020_instant_KmKe_2d_hslice_nt_', string(nt), ...
        '_vlevel_', string(vlevel), '.nc');

elseif(strcmp(case_type, 'avg'))
    % read files and calculate mean (bar) values of u, v
    % indxRange = 469:3382; % entire year
    % indxRange = 712:768; % April 01 - April 07, 2019
    % indxRange = 1040:1096; % May 12 - May 19, 2019
    % indxRange = 2104:2160; % Sept 22 - Sept 29, 2019
    indxRange = 3016:3072; % Jan 14 - Jan 20, 2020

    nt0=indxRange(1);
    [~, Nt] = size(indxRange);
    time_arr = zeros(Nt, 1);
    % ub=zeros(N, Mu, Lu);
    % vb=zeros(N, Mv, Lv);

    KmKe_2d_indxRange = indxRange;
    [~, mean_Nt] = size(KmKe_2d_indxRange);
    KmKe_2d_file = strcat('nesm_2019_2020_KmKe_2d_hslice_nt_', string(KmKe_2d_indxRange(1)), '_', string(KmKe_2d_indxRange(mean_Nt)), ...
        '_vlevel_', string(vlevel), '.nc');
end

HRS2D = ncread(KmKe_2d_file, 'HRS2D');
VRS2D = ncread(KmKe_2d_file, 'VRS2D');
KmKe2D = ncread(KmKe_2d_file, 'KmKe2D');
%--------------------------------------------------------------------------
% only to get the dates in the title!
vort_vlevel = -4500;
vort_indxRange = 0:3382; 
[~, vort_Nt] = size(vort_indxRange);
nt0=vort_indxRange(1);

vort_t_filename = strcat('../nesm_2019_2020_vort_t_hslice_vlevel_', string(vort_vlevel(1)), '_nt_',string(vort_indxRange(1)), '_', ...
    string(vort_indxRange(vort_Nt)), '.nc');
t_arr = ncread(vort_t_filename, 't_arr');
t_i = t_arr(indxRange(1) - nt0);
t_e = t_arr(indxRange(Nt) - nt0);
date_i = char(string(datetime(2017, 12, 31, 23, 29, t_i)));
date_e = char(string(datetime(2017, 12, 31, 23, 29, t_e)));

date_i = date_i(1:6);
date_e = date_e(1:11);

date = strcat(date_i, '-', date_e);
%-------------------------------------------------------

% define box idx to take horizontal averages
box='N' % or 'NE';
[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box);
%--------------------------------------------------------------------------

%% plot the time average, real KmKe2d, now already done in the modified save_KmKe_2d.m script
% KmKe2d_avg = zeros(Ny, Nx);
% for nt=1:Nt
%     KmKe2d_avg = KmKe2d_avg + squeeze(KmKe_t(nt, :, :));
% end
% KmKe2d_avg = KmKe2d_avg./Nt;
%--------------------------------------------------------------------------
%% plot HRS hslice
% figure1 = figure(1);
% % [latlim, lonlim] = geoquadline(lat_rho, lon_rho); geoquadline requires
% % Mapping toolbox
% % Create axes
% ax1 = axes('Parent', figure1, 'YMinorTick','on',...
%     'LineWidth',3,...
%     'FontSize',24);
% 
% m_proj('miller','long', lonlim,'lat', latlim);
% 
% 
% hold on;
% zMin = -5e-7; % min(min(vort)); 
% zMax = -zMin;
% h1 = m_image(lon_rho_vec, lat_rho_vec, HRS2D);
% colormap('parula'); 
% cb = colorbar; 
% cb.FontSize = 20;
% clim([zMin zMax]);
% % colormap(ax1,b2r(zMin,zMax));  colorbar;
% % colormap(ax1, "jet"); clim([zMin zMax]); colorbar;
% % colormap(ax1, whitejet); clim([zMin zMax]); colorbar;
% % freezeColors; hold on;
% 
% [h2, c2] = m_contourf(lon_rho, lat_rho, -depth, [vlevel vlevel], ...
%     'LineWidth', 3, 'EdgeColor', [0 0 0], 'FaceColor', [1  1 1], 'FaceAlpha', 1);
% % c2.FaceColor = [1 1 1]; c2.FaceAlpha = 0.3; % opacity
% % c2.EdgeColor = [0 0 0];
% 
% m_grid('tickdir','in', ...
% 'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
% 'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
% 'ytick',([38 39]), ... % latitude        
% 'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;
% 
% % title(ax1, string(datetime(2017,12,31,0,0,t_arr(vort_t_nt))));
% % ax1.TitleHorizontalAlignment = 'center'; % left makes it come to center
% 
% %%%
% set(figure1, 'Visible', 'off'); % stop pop-ups
% if(strcmp(case_type, 'instant'))
%     figname  = [plots_path, 'rs_plots/hrs/nesm_2019_2020_instant_hrs_2d_hslice'];
% 
%     figname = strcat(figname, '_vlevel_', string(vlevel));
%     figname = strcat(figname, '_nt_', string(nt), ...
%         '_vlevel_', string(vlevel));
% elseif(strcmp(case_type, 'avg'))
%     figname  = [plots_path, 'rs_plots/hrs/nesm_2019_2020_hrs_2d_hslice'];
% 
%     figname = strcat(figname, '_vlevel_', string(vlevel));
%     figname = strcat(figname, '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), ...
%         '_vlevel_', string(vlevel));
% end
% % vort_contour = strcat(vort_contour, '.pdf');
% % exportgraphics(figure1, vort_contour, 'ContentType', 'vector'); % remove extra white space, 2022a and above
% % 
% exportgraphics(figure1,strcat(figname, '.eps'))
% close all;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% plot VRS hslice
% figure1 = figure(2);
% % [latlim, lonlim] = geoquadline(lat_rho, lon_rho); geoquadline requires
% % Mapping toolbox
% % Create axes
% ax1 = axes('Parent', figure1, 'YMinorTick','on',...
%     'LineWidth',3,...
%     'FontSize',24);
% 
% m_proj('miller','long', lonlim,'lat', latlim);
% 
% 
% hold on;
% zMin = -5e-7; % min(min(vort)); 
% zMax = -zMin;
% h1 = m_image(lon_rho_vec, lat_rho_vec, KmKe2D);
% colormap('parula'); 
% cb = colorbar; 
% cb.FontSize = 20;
% clim([zMin zMax]);
% % colormap(ax1,b2r(zMin,zMax));  colorbar;
% % colormap(ax1, "jet"); clim([zMin zMax]); colorbar;
% % colormap(ax1, whitejet); clim([zMin zMax]); colorbar;
% % freezeColors; hold on;
% 
% [h2, c2] = m_contourf(lon_rho, lat_rho, -depth, [vlevel vlevel], ...
%     'LineWidth', 3, 'EdgeColor', [0 0 0], 'FaceColor', [1  1 1], 'FaceAlpha', 1);
% % c2.FaceColor = [1 1 1]; c2.FaceAlpha = 0.3; % opacity
% % c2.EdgeColor = [0 0 0];
% 
% m_grid('tickdir','in', ...
% 'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
% 'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
% 'ytick',([38 39]), ... % latitude        
% 'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;
% 
% % title(ax1, string(datetime(2017,12,31,0,0,t_arr(vort_t_nt))));
% % ax1.TitleHorizontalAlignment = 'center'; % left makes it come to center
% %%%
% set(figure1, 'Visible', 'off'); % stop pop-ups
% if(strcmp(case_type, 'instant'))
%     figname  = [plots_path, 'rs_plots/vrs/nesm_2019_2020_instant_vrs_2d_hslice'];
% 
%     figname = strcat(figname, '_vlevel_', string(vlevel));
%     figname = strcat(figname, '_nt_', string(nt), ...
%         '_vlevel_', string(vlevel));
% elseif(strcmp(case_type, 'avg'))
%     figname  = [plots_path, 'rs_plots/vrs/nesm_2019_2020_vrs_2d_hslice'];
% 
%     figname = strcat(figname, '_vlevel_', string(vlevel));
%     figname = strcat(figname, '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), ...
%         '_vlevel_', string(vlevel));
% end% vort_contour = strcat(vort_contour, '.pdf');
% % exportgraphics(figure1, vort_contour, 'ContentType', 'vector'); % remove extra white space, 2022a and above
% % 
% exportgraphics(figure1,strcat(figname, '.eps'))
% close all;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% plot KmKe hslice
% figure1 = figure(3);
% % [latlim, lonlim] = geoquadline(lat_rho, lon_rho); geoquadline requires
% % Mapping toolbox
% % Create axes
% ax1 = axes('Parent', figure1, 'YMinorTick','on',...
%     'LineWidth',3,...
%     'FontSize',24);
% 
% m_proj('miller','long', lonlim,'lat', latlim);
% 
% 
% hold on;
% zMin = -5e-7; % min(min(vort)); 
% zMax = -zMin;
% h1 = m_image(lon_rho_vec, lat_rho_vec, KmKe2D);
% colormap(cmocean('haline')); 
% colorbar; 
% clim([zMin zMax]);
% % colormap(ax1,b2r(zMin,zMax));  colorbar;
% % colormap(ax1, "jet"); clim([zMin zMax]); colorbar;
% % colormap(ax1, whitejet); clim([zMin zMax]); colorbar;
% % freezeColors; hold on;
% 
% [h2, c2] = m_contourf(lon_rho, lat_rho, -depth, [vlevel vlevel], ...
%     'LineWidth', 3, 'EdgeColor', [0 0 0], 'FaceColor', [1  1 1], 'FaceAlpha', 1);
% % c2.FaceColor = [1 1 1]; c2.FaceAlpha = 0.3; % opacity
% % c2.EdgeColor = [0 0 0];
% 
% m_grid('tickdir','in', ...
% 'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
% 'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
% 'ytick',([38 39]), ... % latitude        
% 'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;
% 
% % title(ax1, string(datetime(2017,12,31,0,0,t_arr(vort_t_nt))));
% % ax1.TitleHorizontalAlignment = 'center'; % left makes it come to center
% hold on;
% % plot boxes over which horizontal averages will be taken
% m_plot([lon_rho_vec(lon_idx_min) lon_rho_vec(lon_idx_max)],[lat_rho_vec(lat_idx_min) lat_rho_vec(lat_idx_min)], ...
%     'r', 'LineWidth', 3);
% m_plot([lon_rho_vec(lon_idx_min) lon_rho_vec(lon_idx_max)],[lat_rho_vec(lat_idx_max) lat_rho_vec(lat_idx_max)], ...
%     'r', 'LineWidth', 3);
% 
% m_plot([lon_rho_vec(lon_idx_min) lon_rho_vec(lon_idx_min)],[lat_rho_vec(lat_idx_min) lat_rho_vec(lat_idx_max)], ...
%     'r', 'LineWidth', 3);
% m_plot([lon_rho_vec(lon_idx_max) lon_rho_vec(lon_idx_max)],[lat_rho_vec(lat_idx_min) lat_rho_vec(lat_idx_max)], ...
%     'r', 'LineWidth', 3);
% %%%
% set(figure1, 'Visible', 'off'); % stop pop-ups
% figname  = [plots_path, 'rs_plots/KmKe/nesm_2019_2020_KmKe_2d_hslice'];
% 
% figname = strcat(figname, '_vlevel_', string(vlevel));
% figname = strcat(figname, '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), ...
%     '_vlevel_', string(vlevel));
% % vort_contour = strcat(vort_contour, '.pdf');
% % exportgraphics(figure1, vort_contour, 'ContentType', 'vector'); % remove extra white space, 2022a and above
% % 
% exportgraphics(figure1,strcat(figname, '.eps'))
% close all;
%--------------------------------------------------------------------------
%% plot KmKe hslice with all the boxes
% box_arr = ["N", "NE", "W", "SW", "S", "SE", "E", "NW"]
% box_arr = ["big_N", "big_S"]
box_arr = []
figure1 = figure(3);
% [latlim, lonlim] = geoquadline(lat_rho, lon_rho); geoquadline requires
% Mapping toolbox
% Create axes
ax1 = axes('Parent', figure1, 'YMinorTick','on',...
    'LineWidth',3,...
    'FontSize',24);

m_proj('miller','long', lonlim,'lat', latlim);


hold on;
zMin = -5e-7; % min(min(vort)); 
zMax = -zMin;
h1 = m_image(lon_rho_vec, lat_rho_vec, KmKe2D);
colormap('parula'); 
cb = colorbar; 
cb.FontSize = 20;
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

title(ax1, string(date));
% title(ax1, string(datetime(2017,12,31,0,0,t_arr(vort_t_nt))));
% ax1.TitleHorizontalAlignment = 'center'; % left makes it come to center
hold on;
% plot boxes over which horizontal averages will be taken
for i = 1:length(box_arr)
    box = string(box_arr(i))
    [lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box);
    
m_plot([lon_rho_vec(lon_idx_min) lon_rho_vec(lon_idx_max)],[lat_rho_vec(lat_idx_min) lat_rho_vec(lat_idx_min)], ...
    'r', 'LineWidth', 3);
m_plot([lon_rho_vec(lon_idx_min) lon_rho_vec(lon_idx_max)],[lat_rho_vec(lat_idx_max) lat_rho_vec(lat_idx_max)], ...
    'r', 'LineWidth', 3);

m_plot([lon_rho_vec(lon_idx_min) lon_rho_vec(lon_idx_min)],[lat_rho_vec(lat_idx_min) lat_rho_vec(lat_idx_max)], ...
    'r', 'LineWidth', 3);
m_plot([lon_rho_vec(lon_idx_max) lon_rho_vec(lon_idx_max)],[lat_rho_vec(lat_idx_min) lat_rho_vec(lat_idx_max)], ...
    'r', 'LineWidth', 3);

hold on;
end
%%%
set(figure1, 'Visible', 'off'); % stop pop-ups
if(strcmp(case_type, 'instant'))
    figname  = [plots_path, 'rs_plots/KmKe/nesm_2019_2020_instant_KmKe_2d_hslice_boxes'];

    figname = strcat(figname, '_vlevel_', string(vlevel));
    figname = strcat(figname, '_nt_', string(nt), ...
    '_vlevel_', string(vlevel));

elseif(strcmp(case_type, 'avg'))
    figname  = [plots_path, 'rs_plots/KmKe/nesm_2019_2020_KmKe_2d_hslice_boxes'];

    figname = strcat(figname, '_vlevel_', string(vlevel));
    figname = strcat(figname, '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), ...
    '_vlevel_', string(vlevel));
end

exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
% export_fig(figure1, figname, '-eps','-transparent', '-r300'); % 
close all;
%--------------------------------------------------------------------------