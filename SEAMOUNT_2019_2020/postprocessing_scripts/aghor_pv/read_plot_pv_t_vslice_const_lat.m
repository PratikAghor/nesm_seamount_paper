clear all; clc;
%-------------------------------------------------------
%%
start_paths
%-------------------------------------------------------
Omega = 7.2921e-5; % rotation rate of Earth s^-1
N = 1e-3; % Brunt-Vaisala frequency (order of magnitude, s^-1)
%%Loop to calculate averages if need (modify as needed)
indxRange = 3016:3072; % Jan 14 - Jan 20
% indxRange = 3016:3024; % Jan 14 - Jan 20
% indxRange = 3024:3032; % 
% indxRange = 3032:3040; % 
% indxRange = 3040:3048; % 
% indxRange = 3048:3056; % 
% indxRange = 3056:3064; % 
% indxRange = 3064:3072; % 
nt0=indxRange(1);
% Aghor's script
% get avg ke as a function of height and temperature
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);

rho_indxRange = 1691:3382;
rho_nt0 = rho_indxRange(1);
[~, rho_Nt] = size(rho_indxRange);

tskip = 1;
vlevel_skip = 50;
%-------------------------------------------------------
%-------------------------------------------------------
%%
pv_t_filename = strcat('seamount_2019_2020_pv_t_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
% rho_pot_t_filename = strcat('seamount_2019_2020_rho_pot_t_vslice_const_lat_nt_',string(rho_indxRange(1)), '_', ...
%     string(rho_indxRange(rho_Nt)),'.nc');

% t_arr = ncread(vort_t_filename, 't_arr');
lat_idx_arr = ncread(pv_t_filename, 'lat_idx');
lat_arr = ncread(pv_t_filename, 'lat');
pv_t_vslice = ncread(pv_t_filename, 'pv_t_vslice');
% rho_pot_t_vslice = ncread(rho_pot_t_filename, 'rho_t_vslice');

sec = 'nozoom' % zoom or nozoom
%%
%-------------------------------------------------------
% only to get the dates in the title!
vlevel = -4500;
vort_indxRange = 0:3382; 
[~, vort_Nt] = size(vort_indxRange);
nt0=vort_indxRange(1);

vort_t_filename = strcat('../seamount_2019_2020_vort_t_hslice_vlevel_', string(vlevel(1)), '_nt_',string(vort_indxRange(1)), '_', ...
    string(vort_indxRange(vort_Nt)), '.nc');
t_arr = ncread(vort_t_filename, 't_arr');
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

% mkdir vort_plots;
%%
% 
%% plotting
for nt = 1:Nt
    vort_t_nt = indxRange(nt) - nt0;
    for i = 1:length(lat_arr)
        % lat_idx = find(abs(lat_rho_vec-lat_arr(1))<1e-3); % find idx of lat in lat_rho_vec
    
        lat_idx = round(lat_idx_arr(i));
        figure1 = figure(i);
        % [latlim, lonlim] = geoquadline(lat_rho, lon_rho); % geoquadline requires Mapping Toolbox.
    
        X = (repmat(lon_rho_vec, NumLayers, 1));
        X = (X);
        Z = (squeeze(zr(:, lat_idx, :)));
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        hold on;
        epv_vslice = squeeze(pv_t_vslice(nt, :, i, :));
        zMin = -3e-9; % min(min(epv)); 
        zMax = -zMin; % max(max(~isinf(w_vslice)));
        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
        % pcolor(ax1, Y(2:end-1, 2:end-1), Z(2:end-1, 2:end-1), epv_vslice(2:end-1, 2:end-1));
        pcolor(ax1, X, Z, epv_vslice); 
        shading interp;
        % freezeColors; hold on;
        set(ax1,'Color', [1 1 1])
        % set(ax1,'YDir','reverse')
        % colormap(ax1,b2r(zMin,zMax));  
        % colormap(ax1, whitejet); clim([zMin zMax]);
        % colormap(ax1, "jet"); 
        colormap(cmocean('balance')); colorbar; clim([zMin zMax]);
        % clim([zMin zMax]);
        colorbar;
        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        title(ax1, date);
        % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        
        hold on;
        % Z = -flip(squeeze(zr(lat_idx,:,:)), 1);
        % LevelList1 = linspace(min(min(rho_vslice)), 36.49, 10);
        % LevelList2 = linspace(36.49, 36.8, 5);
        % LevelList3 = linspace(36.8, max(max(rho_vslice)), 10);
        % 
        % LevelList = [LevelList1, LevelList2, LevelList3];
        % contour(ax1, X, (Z), (rho_vslice), ...
        %     'EdgeColor', [0 0 0], 'LevelList', LevelList);
    
    
        % visual check for Symmetric instability (SI)
        if (strcmp(sec, 'zoom')==true)
            xlim([-63.5, -62.9]);
            ylim([-5000, -4000]);
        elseif(strcmp(sec, 'nozoom')==true)
            xlim([-64., -62.]);
            ylim([-5000, -1500]);
        end
        xticks([-64.99 -64 -63 -62 -61])  % longitude   
        xticklabels({'65°W', '64°W','63°W','62°W','61°W'}) % name longitude ticks as you want
        yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
        yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
        axis square
        % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        %%%
        set(figure1, 'Visible', 'off'); % stop pop-ups
        if (strcmp(sec, 'zoom')==true)
            figname  = [plots_path, 'si_check/vslice/const_lat_zoom/seamount_2019_2020_pv_t_vslice'];
        elseif(strcmp(sec, 'nozoom')==true)
            figname  = [plots_path, 'si_check/vslice/const_lat/seamount_2019_2020_pv_t_vslice'];
        end
        % 
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(indxRange(nt)));
        % vort_contour = strcat(vort_contour, '.pdf');
        % exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector'); % remove extra white space, 2022a and above
        % 
        % set(figure1, 'PaperPositionMode', 'auto')
        % print(figure1,strcat(figname, '.png'),'-dpng','-r300');
        % exportgraphics(figure1, strcat(figname, '.eps'), 'Resolution',300); % remove extra white space, 2022a and above
        export_fig(ax1, figname, '-eps','-transparent', '-r300'); % 
        close all;
    end
end
%------------------------------------

