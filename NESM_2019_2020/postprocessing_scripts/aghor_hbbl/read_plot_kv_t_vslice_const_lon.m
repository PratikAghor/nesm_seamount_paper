clear all; clc;
%-------------------------------------------------------
%%
start_paths
%-------------------------------------------------------
Omega = 7.2921e-5; % rotation rate of Earth s^-1
N = 1e-3; % Brunt-Vaisala frequency (order of magnitude, s^-1)
%%Loop to calculate averages if need (modify as needed)
% indxRange = 469:3382; % Jan 14 - Jan 20
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
% kv_t_filename = strcat('nesm_2019_2020_kv_t_vslice_const_lon_nt_',string(indxRange(1)), '_', ...
%     string(indxRange(Nt)),'_oka_niwa.nc');
kv_t_filename = strcat('nesm_2019_2020_kv_t_vslice_const_lon_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'_laurent_ferrari.nc');
% rho_pot_t_filename = strcat('nesm_2019_2020_rho_pot_t_vslice_const_lat_nt_',string(rho_indxRange(1)), '_', ...
%     string(rho_indxRange(rho_Nt)),'.nc');

% t_arr = ncread(vort_t_filename, 't_arr');
lon_idx_arr = ncread(kv_t_filename, 'lon_idx');
lon_arr = ncread(kv_t_filename, 'lon');
kv_t_vslice = ncread(kv_t_filename, 'kv_t_vslice');
% rho_pot_t_vslice = ncread(rho_pot_t_filename, 'rho_t_vslice');

sec = 'nozoom' % zoom or nozoom
%%
%-------------------------------------------------------
rho_pot_t_filename = strcat('../aghor_pv/nesm_2019_2020_rho_pot_t_vslice_const_lon_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
% rho_pot_t_filename = strcat('nesm_2019_2020_rho_pot_t_vslice_const_lat_nt_',string(rho_indxRange(1)), '_', ...
%     string(rho_indxRange(rho_Nt)),'.nc');

rho_pot_t_vslice = ncread(rho_pot_t_filename, 'rho_pot_t_vslice');
% rho_pot_t_vslice = ncread(rho_pot_t_filename, 'rho_t_vslice');

%%
%-------------------------------------------------------

% only to get the dates in the title!
t_arr = ncread(kv_t_filename, 't_arr');
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
for nt = 28:30
    vort_t_nt = indxRange(nt) - nt0 + 1;
    for j = 1:length(lon_arr)
        % lon_idx = find(abs(lon_rho_vec-lon_arr(j))<1e-3); % find idx of lat in lon_rho_vec
        lon_idx = round(lon_idx_arr(j))        
        
        figure1 = figure(j);
        % [latlim, lonlim] = geoquadline(lat_rho, lon_rho); % geoquadline requires Mapping Toolbox.
    
        Y = repmat(lat_rho_vec', NumLayers, 1);
        Z = (squeeze(zr(:, :, lon_idx)));
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        hold on;
        kv_vslice = squeeze(kv_t_vslice(nt, :, :, j));
        rho_pot_vslice = squeeze(rho_pot_t_vslice(nt, :, :, j));

        zMin = 1e-5; % min(min(epv)); 
        zMax = 1e-2; % max(max(~isinf(w_vslice)));
        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
        % pcolor(ax1, Y(2:end-1, 2:end-1), Z(2:end-1, 2:end-1), epv_vslice(2:end-1, 2:end-1));
        pcolor(ax1, Y, Z, kv_vslice); 
        shading interp;
        % freezeColors; hold on;
        set(ax1,'Color', [1 1 1])
        % colormap(cmocean('balance')); colorbar; clim([zMin zMax]);
        % colormap('sky'); cb = colorbar;
        colormap('jet'); cb = colorbar;
        cb.FontSize = 20; 
        set(ax1,'ColorScale','log')
        clim([zMin zMax]);
        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        % title(ax1, date);
        % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        
        hold on;
        LevelList1 = linspace(min(min(rho_pot_vslice)), 36.49, 20);
        LevelList2 = linspace(36.49, 36.8, 10);
        LevelList3 = linspace(36.8, max(max(rho_pot_vslice)), 20);

        LevelList = [LevelList1, LevelList2, LevelList3];
        contour(ax1, Y, (Z), (rho_pot_vslice), ...
            'EdgeColor', [1 1 1], 'LevelList', LevelList);
    
    
        % visual check for Symmetric instability (SI)
        if (strcmp(sec, 'zoom')==true)
            xlim([38, 39]);
            ylim([-5000, -4000]);
        elseif(strcmp(sec, 'nozoom')==true)
            xlim([latlim]);
            ylim([-5000, -1500]);
        end
        xticks([38 39 39.99])  % longitude   
        xticklabels({'38°N', '39°N','40°N'}) % name longitude ticks as you want
        yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
        yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
        axis square
        % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        %%%
        set(figure1, 'Visible', 'off'); % stop pop-ups
        if (strcmp(sec, 'zoom')==true)
            figname  = [plots_path, 'kv_plots/vslice/const_lon_zoom/nesm_2019_2020_kv_t_vslice'];
        elseif(strcmp(sec, 'nozoom')==true)
            figname  = [plots_path, 'kv_plots/vslice/const_lon/nesm_2019_2020_kv_t_vslice'];
        end
        % 
        figname = strcat(figname, '_lon_', string(lon_arr(j)));
        figname = strcat(figname, '_nt_', string(indxRange(nt)));
        exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector'); % remove extra white space, 2022a and above
        % 
        % set(figure1, 'PaperPositionMode', 'auto')
        % print(figure1,strcat(figname, '.png'),'-dpng','-r300');
        % exportgraphics(figure1, strcat(figname, '.eps'), 'Resolution',300); % remove extra white space, 2022a and above
        % export_fig(ax1, figname, '-eps','-transparent', '-r300'); % 
        close all;
    end
end
%------------------------------------

